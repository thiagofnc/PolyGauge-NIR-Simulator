"""Tests for the shared pipeline, loaders, Excel export and UI/core/export consistency."""

import os
import tempfile
import unittest

import numpy as np

from BroadbandExcel import export_broadband_workbook
from BroadbandPipeline import (DETECTOR_RESPONSIVITY, LAYERS_PHYSICAL, LAYERS_RATIO, REFERENCE_MULTIPLIER,
                               SOURCE_BLACKBODY, BroadbandConfig, run_broadband_analysis)
from SpectralData import load_material_library, load_source_presets, load_spectral_curve_csv, load_material_spectrum

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def synthetic_library(absorbance=0.2, low=2000.0, high=3000.0, reference_um=None):
    wl = np.linspace(low, high, 201)
    values = absorbance(wl) if callable(absorbance) else np.full_like(wl, absorbance)
    return {"Film": {"name": "Film", "wavelength_nm": wl, "absorbance": values, "convention": "base10",
                     "reference_thickness_um": reference_um, "source": "synthetic", "note": ""}}


def config(**overrides):
    base = dict(material="Film", values=[1, 2, 3, 4], no_film_voltage_mv=100.0, layer_to_reference_ratio=1.0,
                wavelength_min_nm=2000.0, wavelength_max_nm=3000.0)
    base.update(overrides)
    return BroadbandConfig(**base)


class PipelineTests(unittest.TestCase):
    def test_corrections_are_off_by_default(self):
        cfg = config()
        self.assertFalse(cfg.baseline_correction or cfg.clamp_negative or cfg.interface_correction)
        self.assertEqual(cfg.dark_voltage_mv, 0.0)
        result = run_broadband_analysis(cfg, synthetic_library())
        self.assertFalse(result["corrections"]["any_enabled"])
        for row in result["rows"]:
            self.assertEqual(row["predicted_voltage_mv"], row["raw_predicted_voltage_mv"])

    def test_physical_thickness_uses_x_over_x_ref(self):
        result = run_broadband_analysis(config(thickness_mode=LAYERS_PHYSICAL, layer_thickness_um=25.0),
                                        synthetic_library(0.2, reference_um=50.0))
        row = result["rows"][3]  # 4 layers x 25 µm = 100 µm = 2 x_ref
        self.assertEqual(row["thickness_um"], 100.0)
        self.assertEqual(row["thickness_scale"], 2.0)
        self.assertAlmostEqual(row["predicted_transmission"], 10 ** (-2 * 0.2), places=12)
        self.assertTrue(result["thickness"]["reports_um"])

    def test_layer_thickness_value_changes_prediction(self):
        library = synthetic_library(0.2, reference_um=50.0)
        thin = run_broadband_analysis(config(thickness_mode=LAYERS_PHYSICAL, layer_thickness_um=10.0), library)
        thick = run_broadband_analysis(config(thickness_mode=LAYERS_PHYSICAL, layer_thickness_um=40.0), library)
        self.assertGreater(thin["rows"][0]["predicted_voltage_mv"], thick["rows"][0]["predicted_voltage_mv"])

    def test_physical_mode_refuses_unknown_reference_thickness(self):
        with self.assertRaises(ValueError):
            run_broadband_analysis(config(thickness_mode=LAYERS_PHYSICAL, layer_thickness_um=25.0), synthetic_library())

    def test_ratio_mode_never_reports_micrometres(self):
        result = run_broadband_analysis(config(layer_to_reference_ratio=0.5), synthetic_library(0.2))
        self.assertFalse(result["thickness"]["reports_um"])
        self.assertIsNone(result["rows"][0]["thickness_um"])
        self.assertAlmostEqual(result["rows"][3]["predicted_transmission"], 10 ** (-0.4), places=12)

    def test_dark_voltage_in_pipeline(self):
        result = run_broadband_analysis(config(dark_voltage_mv=12.0, measured_voltages_by_layer={1: 60.0}),
                                        synthetic_library(0.2))
        row = result["rows"][0]
        self.assertAlmostEqual(row["predicted_voltage_mv"], 12.0 + 88.0 * 10 ** -0.2, places=12)
        self.assertAlmostEqual(row["measured_transmission"], (60.0 - 12.0) / 88.0, places=12)

    def test_measured_layers_not_matched_in_multiplier_mode(self):
        result = run_broadband_analysis(config(thickness_mode=REFERENCE_MULTIPLIER,
                                               measured_voltages_by_layer={1: 50.0, 2: 40.0}), synthetic_library())
        self.assertTrue(all(row["measured_voltage_mv"] is None for row in result["rows"]))
        self.assertEqual(result["metrics"]["n"], 0)
        self.assertTrue(any("NOT compared" in w for w in result["warnings"]))

    def test_measured_comparison_and_metrics(self):
        library = synthetic_library(0.2)
        result = run_broadband_analysis(config(values=[1, 2], measured_voltages_by_layer={1: 60.0, 2: 40.0}), library)
        predicted = [10 ** -0.2 * 100, 10 ** -0.4 * 100]
        errors = [predicted[0] - 60.0, predicted[1] - 40.0]
        self.assertAlmostEqual(result["rows"][0]["error_mv"], errors[0], places=10)
        self.assertAlmostEqual(result["rows"][1]["percent_error"], 100 * errors[1] / 40.0, places=10)
        self.assertAlmostEqual(result["metrics"]["mae_mv"], np.mean(np.abs(errors)), places=10)
        self.assertAlmostEqual(result["metrics"]["rmse_mv"], np.sqrt(np.mean(np.square(errors))), places=10)

    def test_inspected_value_is_never_snapped(self):
        with self.assertRaises(ValueError):
            run_broadband_analysis(config(inspect_value=2.2), synthetic_library())
        result = run_broadband_analysis(config(inspect_value=3), synthetic_library())
        self.assertEqual(result["inspected_index"], 2)

    def test_partial_coverage_flagged_and_reported(self):
        result = run_broadband_analysis(config(wavelength_max_nm=4000.0), synthetic_library(0.2))
        self.assertTrue(result["coverage"]["partial"])
        self.assertAlmostEqual(result["coverage"]["fraction"], 0.5, places=12)
        self.assertTrue(result["warnings"][0].startswith("PARTIAL SPECTRAL COVERAGE"))
        low, high = result["rows"][0]["predicted_voltage_bounds_mv"]
        self.assertLess(low, result["rows"][0]["predicted_voltage_mv"])
        self.assertGreater(high, result["rows"][0]["predicted_voltage_mv"])

    def test_negative_absorbance_warning_only_counts_active_region(self):
        library = synthetic_library(lambda wl: np.where(wl < 2500, -0.05, 0.1))
        inside = run_broadband_analysis(config(), library)
        outside = run_broadband_analysis(config(wavelength_min_nm=2600.0), library)
        self.assertGreater(inside["negative_absorbance"]["points"], 0)
        self.assertEqual(outside["negative_absorbance"]["points"], 0)

    def test_enabled_corrections_are_reported(self):
        library = synthetic_library(lambda wl: np.where(wl < 2500, -0.05, 0.1))
        result = run_broadband_analysis(config(clamp_negative=True, interface_correction=True, refractive_index=1.5),
                                        library)
        corrections = result["corrections"]
        self.assertTrue(corrections["any_enabled"])
        self.assertGreater(corrections["clamped_points"], 0)
        self.assertAlmostEqual(corrections["interface_per_layer"], 0.96 ** 2)
        row = result["rows"][1]
        self.assertAlmostEqual(row["interface_transmission"], 0.96 ** 4, places=12)
        self.assertNotEqual(row["predicted_voltage_mv"], row["raw_predicted_voltage_mv"])

    def test_interface_correction_needs_layer_mode(self):
        with self.assertRaises(ValueError):
            run_broadband_analysis(config(thickness_mode=REFERENCE_MULTIPLIER, interface_correction=True,
                                          refractive_index=1.5), synthetic_library())

    def test_responsivity_curve_zero_outside_range(self):
        response = {"wavelength_nm": np.array([2200.0, 2800.0]), "values": np.array([1.0, 1.0]), "note": "test"}
        library = synthetic_library(lambda wl: np.where(wl < 2200, 5.0, 0.1))
        result = run_broadband_analysis(config(detector_mode=DETECTOR_RESPONSIVITY, detector_responsivity=response,
                                               wavelength_min_nm=None, wavelength_max_nm=None), library)
        self.assertAlmostEqual(result["rows"][0]["spectral_transmission"], 10 ** -0.1, places=10)
        self.assertEqual(result["weighting"]["weight_domain_nm"], (2200.0, 2800.0))

    def test_blackbody_source_changes_prediction(self):
        library = synthetic_library(lambda wl: (wl - 2000) / 1000)
        flat = run_broadband_analysis(config(), library)
        hot = run_broadband_analysis(config(source_mode=SOURCE_BLACKBODY, source_temperature_k=903.15), library)
        self.assertNotAlmostEqual(flat["rows"][0]["spectral_transmission"], hot["rows"][0]["spectral_transmission"], 4)
        with self.assertRaises(ValueError):
            run_broadband_analysis(config(source_mode=SOURCE_BLACKBODY), library)


class LoaderTests(unittest.TestCase):
    def test_evoh_loader_converts_fraction_transmission(self):
        entry = {"file": "EVOHTransmissionPercentvsWavelength_in_nm_2000nm_to_5000nm.csv",
                 "format": "csv_wavelength_nm", "quantity": "transmittance_fraction", "reference_thickness_um": None}
        material = load_material_spectrum(BASE_DIR, "EVOH", entry)
        import pandas as pd
        frame = pd.read_csv(os.path.join(BASE_DIR, entry["file"]), skipinitialspace=True).sort_values("x")
        np.testing.assert_allclose(material["wavelength_nm"], frame["x"].to_numpy())
        np.testing.assert_allclose(material["absorbance"], -np.log10(frame["y"].to_numpy()))
        self.assertEqual(material["convention"], "base10")
        self.assertIsNone(material["reference_thickness_um"])

    def test_transmittance_units_are_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.csv")
            with open(path, "w") as handle:
                handle.write("x,y\n2000,90\n3000,50\n")
            entry = {"file": "t.csv", "format": "csv_wavelength_nm", "reference_thickness_um": 30}
            with self.assertRaises(ValueError):
                load_material_spectrum(tmp, "bad", dict(entry, quantity="transmittance_fraction"))
            percent = load_material_spectrum(tmp, "ok", dict(entry, quantity="transmittance_percent"))
            np.testing.assert_allclose(percent["absorbance"], -np.log10([0.9, 0.5]))
            self.assertEqual(percent["reference_thickness_um"], 30.0)

    def test_bundled_library_has_unknown_reference_thickness(self):
        library = load_material_library(BASE_DIR, include_discovered=False)
        self.assertEqual(set(library), {"Nylon", "PE", "EVOH"})
        self.assertTrue(all(m["reference_thickness_um"] is None for m in library.values()))

    def test_curve_csv_units_and_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "r.csv")
            with open(path, "w") as handle:
                handle.write("wavelength (um),responsivity\n3,0.5\n2,0.2\n")
            curve = load_spectral_curve_csv(path, "responsivity")
            np.testing.assert_allclose(curve["wavelength_nm"], [2000.0, 3000.0])
            np.testing.assert_allclose(curve["values"], [0.2, 0.5])
            with open(path, "w") as handle:
                handle.write("wavelength_nm,responsivity\n2000,-0.1\n3000,0.5\n")
            with self.assertRaises(ValueError):
                load_spectral_curve_csv(path, "responsivity")
            with open(path, "w") as handle:
                handle.write("wavelength_nm,T\n2000,80\n3000,50\n")
            with self.assertRaises(ValueError):
                load_spectral_curve_csv(path, "transmission")

    def test_hpir104_preset_from_component_database(self):
        presets = load_source_presets(BASE_DIR)
        hpir = [p for p in presets.values() if p["name"].startswith("HPIR104")]
        self.assertEqual(len(hpir), 1)
        self.assertAlmostEqual(hpir[0]["temperature_k"], 903.15)


class ExportTests(unittest.TestCase):
    def test_excel_export_matches_pipeline(self):
        from openpyxl import load_workbook
        result = run_broadband_analysis(config(wavelength_max_nm=3500.0, dark_voltage_mv=5.0,
                                               measured_voltages_by_layer={1: 60.0, 3: 30.0}, inspect_value=3),
                                        synthetic_library(0.2))
        with tempfile.TemporaryDirectory() as tmp:
            path = export_broadband_workbook(os.path.join(tmp, "out.xlsx"), result, [])
            wb = load_workbook(path)
            self.assertEqual(wb.sheetnames, ["Summary", "Results", "Spectral Calculation", "Layer Predictions",
                                             "Experimental Data", "Charts"])
            results = wb["Results"]
            headers = [cell.value for cell in results[1]]
            predicted_col = headers.index("Predicted Voltage (mV)") + 1
            for i, row in enumerate(result["rows"], start=2):
                self.assertAlmostEqual(results.cell(row=i, column=predicted_col).value, row["predicted_voltage_mv"], places=9)
            self.assertIn("PARTIAL SPECTRAL COVERAGE", str(wb["Summary"]["A4"].value))
            self.assertEqual(wb["Spectral Calculation"].max_row, len(result["grid"]["wavelength_nm"]) + 1)
            self.assertEqual(len(wb["Charts"]._charts), 7)


def _tk_available():
    try:
        import tkinter
        root = tkinter.Tk()
        root.destroy()
        return True
    except Exception:
        return False


@unittest.skipUnless(_tk_available(), "Tk display not available")
class UiConsistencyTests(unittest.TestCase):
    def test_ui_core_export_agree(self):
        from tkinter import messagebox
        import customtkinter as ctk
        from openpyxl import load_workbook
        import BroadbandUI

        errors = []
        original = messagebox.showerror
        messagebox.showerror = lambda title, message, **kw: errors.append(message)
        root = ctk.CTk()
        root.withdraw()
        try:
            window = BroadbandUI.BroadbandAnalysisWindow(root)
            self.assertFalse(window.baseline_var.get() or window.clamp_var.get() or window.interface_var.get())
            window.material_var.set("Nylon")
            window.detector_var.set("Hamamatsu")
            window._detector_changed()
            self.assertEqual(window.v0_var.get(), "498.0")  # detector change must not overwrite V0
            window.apply_detector_defaults()
            window.vdark_var.set("2")
            window.update_analysis()
            self.assertEqual(errors, [])
            core = run_broadband_analysis(window.build_config(), window.materials)
            for ui_row, core_row in zip(window.result["rows"], core["rows"]):
                self.assertEqual(ui_row["predicted_voltage_mv"], core_row["predicted_voltage_mv"])
            with tempfile.TemporaryDirectory() as tmp:
                path = export_broadband_workbook(os.path.join(tmp, "ui.xlsx"), window.result, window.experimental_rows())
                sheet = load_workbook(path)["Results"]
                headers = [cell.value for cell in sheet[1]]
                col = headers.index("Predicted Voltage (mV)") + 1
                for i, row in enumerate(core["rows"], start=2):
                    self.assertAlmostEqual(sheet.cell(row=i, column=col).value, row["predicted_voltage_mv"], places=9)
            for graph in BroadbandUI.GRAPH_TYPES:
                window.graph_var.set(graph)
                window.draw_graph()
            window.inspect_var.set("2.5")
            window.update_analysis()
            self.assertTrue(any("not one of the simulated values" in e for e in errors))
        finally:
            messagebox.showerror = original
            root.destroy()


if __name__ == "__main__":
    unittest.main()
