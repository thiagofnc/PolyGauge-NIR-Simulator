"""Tests for the broadband physics primitives in BroadbandAnalysis.py."""

import math
import unittest

import numpy as np

from BroadbandAnalysis import (absorbance_to_transmission, agreement_metrics, calculate_weighted_transmission,
                               dark_corrected_voltage, estimate_baseline_offset, fresnel_layer_transmission,
                               measured_transmission, predict_detector_voltage, solve_thickness_scale,
                               thickness_scale)


class AbsorbanceConversionTests(unittest.TestCase):
    def setUp(self):
        self.wl = np.linspace(2000.0, 3000.0, 101)

    def test_zero_absorbance_returns_v0(self):
        result = predict_detector_voltage(self.wl, np.zeros_like(self.wl), 498.0)
        self.assertAlmostEqual(result["effective_transmission"], 1.0)
        self.assertAlmostEqual(result["predicted_voltage_mv"], 498.0)

    def test_constant_base10_absorbance_scales_by_layers(self):
        absorbance = np.full_like(self.wl, 0.1)
        one = predict_detector_voltage(self.wl, absorbance, 1.0, 1)
        two = predict_detector_voltage(self.wl, absorbance, 1.0, 2)
        self.assertAlmostEqual(one["effective_transmission"], 10 ** -0.1, places=12)
        self.assertAlmostEqual(two["effective_transmission"], 10 ** -0.2, places=12)
        self.assertNotAlmostEqual(one["effective_transmission"], math.exp(-0.1), places=3)

    def test_constant_natural_attenuation_is_beer_lambert(self):
        result = predict_detector_voltage(self.wl, np.full_like(self.wl, 0.4), 1.0, 2.5, absorbance_mode="natural")
        self.assertAlmostEqual(result["effective_transmission"], np.exp(-1.0), places=12)

    def test_twice_reference_thickness_squares_transmission(self):
        scale = thickness_scale(100.0, 50.0)
        self.assertEqual(scale, 2.0)
        absorbance = np.full_like(self.wl, 0.17)
        result = predict_detector_voltage(self.wl, absorbance, 1.0, scale)
        self.assertAlmostEqual(result["effective_transmission"], 10 ** (-2 * 0.17), places=12)

    def test_physical_thickness_requires_reference(self):
        with self.assertRaises(ValueError):
            thickness_scale(100.0, None)

    def test_clamp_behavior(self):
        absorbance = np.full_like(self.wl, -0.1)
        unclamped = predict_detector_voltage(self.wl, absorbance, 1.0, 1)
        clamped = predict_detector_voltage(self.wl, absorbance, 1.0, 1, clamp_negative=True)
        self.assertAlmostEqual(unclamped["effective_transmission"], 10 ** 0.1, places=12)
        self.assertAlmostEqual(clamped["effective_transmission"], 1.0, places=12)
        self.assertTrue(np.all(clamped["absorbance"] >= 0))
        np.testing.assert_allclose(absorbance_to_transmission(absorbance, 1.0, clamp_negative=True), 1.0)

    def test_dark_voltage_model_exact(self):
        absorbance = np.full_like(self.wl, 0.3)
        result = predict_detector_voltage(self.wl, absorbance, 100.0, 1, dark_voltage_mv=20.0)
        self.assertAlmostEqual(result["predicted_voltage_mv"], 20.0 + (100.0 - 20.0) * 10 ** -0.3, places=12)
        self.assertAlmostEqual(dark_corrected_voltage(0.5, 100.0, 20.0), 60.0, places=12)
        self.assertAlmostEqual(measured_transmission(60.0, 100.0, 20.0), 0.5, places=12)
        with self.assertRaises(ValueError):
            dark_corrected_voltage(0.5, 10.0, 20.0)


class WeightingTests(unittest.TestCase):
    def setUp(self):
        self.wl = np.linspace(2000.0, 3000.0, 101)

    def test_detector_weighting_changes_wavelength_dependent_result(self):
        absorbance = np.linspace(0.0, 1.0, self.wl.size)
        left = (self.wl, np.linspace(1.0, 0.0, self.wl.size))
        right = (self.wl, np.linspace(0.0, 1.0, self.wl.size))
        result_left = predict_detector_voltage(self.wl, absorbance, 1.0, responsivity=left)
        result_right = predict_detector_voltage(self.wl, absorbance, 1.0, responsivity=right)
        self.assertGreater(result_left["predicted_voltage_mv"], result_right["predicted_voltage_mv"])

    def test_non_constant_interpolation_matches_analytic_integral(self):
        # T(λ) = u/1000 on a 10 nm grid, R(λ) = u on a coarse, offset grid, u = λ - 2000.
        transmission = (self.wl - 2000.0) / 1000.0
        response_x = np.array([2000.0, 2255.0, 2730.0, 3000.0])
        response = (response_x, response_x - 2000.0)
        result = calculate_weighted_transmission(self.wl, transmission, responsivity=response)
        # ∫u·u/1000 du / ∫u du over [0, 1000] = 2/3
        self.assertAlmostEqual(result["effective_transmission"], 2.0 / 3.0, places=4)
        index = int(np.where(result["wavelength_nm"] == 2255.0)[0][0])
        self.assertAlmostEqual(result["transmission"][index], 0.255, places=12)
        self.assertAlmostEqual(result["responsivity"][index], 255.0, places=12)

    def test_source_times_detector_weighting(self):
        source = (np.array([2000.0, 3000.0]), np.array([1.0, 2.0]))
        response = (np.linspace(2000.0, 3000.0, 7), np.linspace(2.0, 1.0, 7))
        transmission = (self.wl - 2000.0) / 1000.0
        result = calculate_weighted_transmission(self.wl, transmission, source=source, responsivity=response)
        grid = result["wavelength_nm"]
        u = (grid - 2000.0) / 1000.0
        np.testing.assert_allclose(result["weight"], (1 + u) * (2 - u), rtol=1e-12)
        # ∫(1+u)(2-u)u du / ∫(1+u)(2-u) du over [0,1] = (7/6 + ... ) computed exactly:
        numerator = 1.0 + 1.0 / 3.0 - 0.25   # ∫(2u + u² - u³) du
        denominator = 2.0 + 0.5 - 1.0 / 3.0  # ∫(2 + u - u²) du
        self.assertAlmostEqual(result["effective_transmission"], numerator / denominator, places=4)

    def test_multiplicative_optical_term_selects_band(self):
        transmission = np.where(self.wl <= 2500.0, 0.2, 0.8)
        filter_curve = (np.array([2000.0, 2500.0, 2500.0001, 3000.0]), np.array([1.0, 1.0, 0.0, 0.0]))
        with_filter = calculate_weighted_transmission(self.wl, transmission, multiplicative_terms=[filter_curve])
        without = calculate_weighted_transmission(self.wl, transmission)
        self.assertAlmostEqual(with_filter["effective_transmission"], 0.2, places=3)
        self.assertAlmostEqual(without["effective_transmission"], 0.5, places=2)
        constant_filter = (self.wl, np.full_like(self.wl, 0.5))
        scaled = calculate_weighted_transmission(self.wl, transmission, multiplicative_terms=[constant_filter])
        self.assertAlmostEqual(scaled["effective_transmission"], without["effective_transmission"], places=12)

    def test_integration_is_restricted_to_overlap(self):
        response = (np.array([2400.0, 2500.0, 2600.0]), np.ones(3))
        result = calculate_weighted_transmission(self.wl, np.ones_like(self.wl), responsivity=response,
                                                 wavelength_min=2000.0, wavelength_max=3000.0)
        self.assertEqual(result["integration_range_nm"], (2400.0, 2600.0))
        self.assertGreaterEqual(result["wavelength_nm"].min(), 2400.0)
        self.assertLessEqual(result["wavelength_nm"].max(), 2600.0)
        self.assertEqual(result["coverage_fraction"], 1.0)

    def test_incomplete_material_coverage_is_flagged(self):
        result = predict_detector_voltage(self.wl, np.full_like(self.wl, 0.3), 100.0, 1,
                                          wavelength_min=2000.0, wavelength_max=4000.0)
        self.assertAlmostEqual(result["coverage_fraction"], 0.5, places=12)
        self.assertTrue(result["partial_coverage"])
        t = 10 ** -0.3
        self.assertAlmostEqual(result["spectral_transmission"], t, places=12)
        low, high = result["effective_transmission_bounds"]
        self.assertAlmostEqual(low, 0.5 * t, places=12)
        self.assertAlmostEqual(high, 0.5 * t + 0.5, places=12)
        self.assertTrue(np.all(np.isnan(result["transmission"][result["wavelength_nm"] > 3000.0])))

    def test_flat_weighting_without_window_uses_material_range(self):
        result = calculate_weighted_transmission(self.wl, np.full_like(self.wl, 0.4))
        self.assertEqual(result["integration_range_nm"], (2000.0, 3000.0))

    def test_bare_arrays_are_rejected_for_weighting_curves(self):
        with self.assertRaises(ValueError):
            calculate_weighted_transmission(self.wl, np.ones_like(self.wl), responsivity=np.ones_like(self.wl))


class CorrectionAndMetricTests(unittest.TestCase):
    def test_agreement_metrics(self):
        metrics = agreement_metrics([1.0, 2.0, 3.0], [1.0, 3.0, 2.0])
        self.assertAlmostEqual(metrics["mae_mv"], 2.0 / 3.0, places=12)
        self.assertAlmostEqual(metrics["rmse_mv"], math.sqrt(2.0 / 3.0), places=12)
        self.assertAlmostEqual(metrics["r_squared"], 0.0, places=12)
        self.assertAlmostEqual(metrics["mape_percent"], 100 * (0 + 0.5 + 1 / 3) / 3, places=10)
        perfect = agreement_metrics([10.0, 20.0], [10.0, 20.0])
        self.assertEqual(perfect["r_squared"], 1.0)
        self.assertLess(agreement_metrics([1.0, 2.0], [3.0, 0.0])["r_squared"], 0.0)

    def test_baseline_offset_ignores_out_of_band_values(self):
        wl = np.linspace(2000.0, 3000.0, 101)
        absorbance = np.full_like(wl, -0.03)
        absorbance[40:45] = 0.4
        absorbance[wl < 2100] = -1.0
        self.assertAlmostEqual(estimate_baseline_offset(wl, absorbance, 2200, 3000), -0.03)

    def test_fresnel_loss(self):
        self.assertAlmostEqual(fresnel_layer_transmission(1.5), 0.96 ** 2)


class ThicknessSolverTests(unittest.TestCase):
    """The inverse model must invert the forward model exactly, not approximately."""

    def setUp(self):
        self.wl = np.linspace(2000.0, 3000.0, 201)
        self.absorbance = np.full_like(self.wl, 0.3)

    def forward(self, scale):
        return predict_detector_voltage(self.wl, self.absorbance, 100.0, scale)["effective_transmission"]

    def test_round_trip_recovers_the_thickness_scale(self):
        for expected in (0.25, 1.0, 2.5, 7.0):
            solution = solve_thickness_scale(self.forward, self.forward(expected))
            self.assertEqual(solution["status"], "ok")
            self.assertTrue(solution["converged"])
            self.assertAlmostEqual(solution["thickness_scale"], expected, places=6)

    def test_structured_spectrum_round_trip(self):
        # A wavelength-dependent spectrum: T_eff is an integral with no closed-form inverse.
        self.absorbance = 0.1 + 0.9 * np.exp(-((self.wl - 2400.0) / 80.0) ** 2)
        solution = solve_thickness_scale(self.forward, self.forward(3.3))
        self.assertAlmostEqual(solution["thickness_scale"], 3.3, places=6)

    def test_transmission_at_or_above_no_film_reports_zero_thickness(self):
        solution = solve_thickness_scale(self.forward, 1.0)
        self.assertEqual(solution["status"], "at_or_above_no_film")
        self.assertEqual(solution["thickness_scale"], 0.0)

    def test_target_below_reachable_range_is_flagged_not_guessed(self):
        solution = solve_thickness_scale(self.forward, 1e-12, max_scale=4.0)
        self.assertEqual(solution["status"], "exceeds_max_scale")
        self.assertEqual(solution["thickness_scale"], 4.0)
        self.assertFalse(solution["converged"])

    def test_non_physical_targets_are_rejected(self):
        for target in (0.0, -0.2, float("nan")):
            with self.assertRaises(ValueError):
                solve_thickness_scale(self.forward, target)

    def test_solution_is_consistent_with_measured_transmission(self):
        voltage = 100.0 * self.forward(2.0)
        solution = solve_thickness_scale(self.forward, measured_transmission(voltage, 100.0, 0.0))
        self.assertAlmostEqual(solution["thickness_scale"], 2.0, places=6)
        self.assertAlmostEqual(dark_corrected_voltage(solution["transmission"], 100.0), voltage, places=6)


if __name__ == "__main__":
    unittest.main()
