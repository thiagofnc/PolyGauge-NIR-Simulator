import unittest

import numpy as np

from BroadbandAnalysis import (analyze_layers, calculate_weighted_transmission,
                               estimate_baseline_offset,
                               fresnel_layer_transmission,
                               predict_detector_voltage)


class BroadbandAnalysisTests(unittest.TestCase):
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

    def test_constant_natural_attenuation_is_beer_lambert(self):
        result = predict_detector_voltage(self.wl, np.full_like(self.wl, 0.4), 1.0,
                                          2.5, absorbance_mode="natural")
        self.assertAlmostEqual(result["effective_transmission"], np.exp(-1.0), places=12)

    def test_detector_weighting_changes_wavelength_dependent_result(self):
        absorbance = np.linspace(0.0, 1.0, self.wl.size)
        left = (self.wl, np.linspace(1.0, 0.0, self.wl.size))
        right = (self.wl, np.linspace(0.0, 1.0, self.wl.size))
        result_left = predict_detector_voltage(self.wl, absorbance, 1.0, responsivity=left)
        result_right = predict_detector_voltage(self.wl, absorbance, 1.0, responsivity=right)
        self.assertGreater(result_left["predicted_voltage_mv"], result_right["predicted_voltage_mv"])

    def test_differently_sampled_spectra_are_interpolated(self):
        transmission = np.full_like(self.wl, 0.5)
        response = (np.array([2000.0, 2250.0, 2750.0, 3000.0]), np.ones(4))
        result = calculate_weighted_transmission(self.wl, transmission, responsivity=response)
        self.assertAlmostEqual(result["effective_transmission"], 0.5, places=12)
        self.assertIn(2250.0, result["wavelength_nm"])

    def test_integration_is_restricted_to_overlap(self):
        response = (np.array([2400.0, 2500.0, 2600.0]), np.ones(3))
        result = calculate_weighted_transmission(self.wl, np.ones_like(self.wl), responsivity=response)
        self.assertEqual(result["integration_range_nm"], (2400.0, 2600.0))
        self.assertGreaterEqual(result["wavelength_nm"].min(), 2400.0)
        self.assertLessEqual(result["wavelength_nm"].max(), 2600.0)

    def test_baseline_offset_ignores_out_of_band_values(self):
        absorbance = np.full_like(self.wl, -0.03)
        absorbance[40:45] = 0.4  # absorption peak
        absorbance[self.wl < 2100] = -1.0  # noisy edge outside the detector band
        offset = estimate_baseline_offset(self.wl, absorbance, 2200, 3000)
        self.assertAlmostEqual(offset, -0.03)

    def test_negative_baseline_no_longer_increases_voltage_with_layers(self):
        absorbance = np.full_like(self.wl, -0.03)
        absorbance[40:60] = 0.4
        corrected = absorbance - estimate_baseline_offset(self.wl, absorbance)
        analysis = analyze_layers(self.wl, corrected, 100.0, [1, 2, 3])
        voltages = [row["predicted_voltage_mv"] for row in analysis["rows"]]
        self.assertTrue(voltages[0] > voltages[1] > voltages[2])

    def test_fresnel_loss_applies_once_per_layer(self):
        interface = fresnel_layer_transmission(1.5)
        self.assertAlmostEqual(interface, 0.96 ** 2)
        result = predict_detector_voltage(self.wl, np.zeros_like(self.wl), 100.0, 3,
                                          layer_interface_transmission=interface)
        self.assertAlmostEqual(result["spectral_transmission"], 1.0)
        self.assertAlmostEqual(result["predicted_voltage_mv"], 100.0 * interface ** 3)

    def test_analyze_layers_reports_errors_against_measurements(self):
        analysis = analyze_layers(self.wl, np.zeros_like(self.wl), 200.0, [1, 2],
                                  {1: 250.0, 2: 200.0})
        first, second = analysis["rows"]
        self.assertAlmostEqual(first["error_mv"], -50.0)
        self.assertAlmostEqual(first["percent_error"], -20.0)
        self.assertAlmostEqual(second["error_mv"], 0.0)
        self.assertAlmostEqual(analysis["metrics"]["mape_percent"], 10.0)


if __name__ == "__main__":
    unittest.main()
