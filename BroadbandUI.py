"""CustomTkinter window for broadband detector analysis and export.

This module only collects settings, calls ``run_broadband_analysis`` and draws
its result.  It contains no spectral physics.
"""

from __future__ import annotations

import math
import os
from tkinter import filedialog, messagebox

import customtkinter as ctk
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from BroadbandExcel import export_broadband_workbook
from BroadbandPipeline import (DETECTOR_FLAT_BAND, DETECTOR_RESPONSIVITY, LAYERS_PHYSICAL, LAYERS_RATIO,
                               REFERENCE_MULTIPLIER, SOURCE_BLACKBODY, SOURCE_FLAT, SOURCE_MEASURED,
                               BroadbandConfig, run_broadband_analysis)
from SpectralData import (load_detector_definitions, load_material_library, load_source_presets,
                          load_spectral_curve_csv)

GRAPH_ABSORBANCE = "Absorbance vs wavelength"
GRAPH_TRANSMISSION = "Spectral transmission vs wavelength"
GRAPH_WEIGHTING = "Detector/source weighting vs wavelength"
GRAPH_CONTRIBUTION = "Weighted transmitted contribution"
GRAPH_VOLTAGE = "Predicted voltage vs thickness"
GRAPH_COMPARE = "Measured vs predicted voltage"
GRAPH_PARITY = "Parity plot (measured vs predicted)"
GRAPH_ERROR = "Prediction error vs thickness"
GRAPH_EFFECTIVE_T = "Effective transmission vs thickness"
GRAPH_ATTENUATION = "Effective attenuation ln(1/T) vs thickness"
GRAPH_TYPES = [GRAPH_ABSORBANCE, GRAPH_TRANSMISSION, GRAPH_WEIGHTING, GRAPH_CONTRIBUTION, GRAPH_VOLTAGE,
               GRAPH_COMPARE, GRAPH_PARITY, GRAPH_ERROR, GRAPH_EFFECTIVE_T, GRAPH_ATTENUATION]

THICKNESS_MODES = {
    "Layers × layer thickness (µm, needs x_ref)": LAYERS_PHYSICAL,
    "Layers × thickness ratio (x_ref unknown)": LAYERS_RATIO,
    "Reference-sample multiplier": REFERENCE_MULTIPLIER,
}
DETECTOR_MODES = {"Flat band approximation": DETECTOR_FLAT_BAND,
                  "Loaded detector responsivity (CSV)": DETECTOR_RESPONSIVITY}
SOURCE_MODES = {"Flat (no source weighting)": SOURCE_FLAT,
                "Blackbody approximation": SOURCE_BLACKBODY,
                "Measured source spectrum (CSV)": SOURCE_MEASURED}
CONVENTIONS = {"From material data": None, "Override: base-10 absorbance": "base10",
               "Override: natural-log attenuation": "natural"}
CUSTOM_TEMPERATURE = "Custom temperature"
TYPICAL_INDEX_NOTE = "Typical mid-IR n: PE ≈ 1.51, EVOH ≈ 1.52, Nylon ≈ 1.53"
SERIES_COLORS = ["#38bdf8", "#34d399", "#fbbf24", "#fb7185", "#a78bfa", "#f472b6", "#2dd4bf", "#fb923c"]


def _optional_float(text, name):
    text = str(text).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{name} must be a number (got {text!r}).") from None


def _required_float(text, name):
    value = _optional_float(text, name)
    if value is None:
        raise ValueError(f"{name} is required.")
    return value


class BroadbandAnalysisWindow(ctk.CTkToplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        self.title("Broadband Spectral Detector Analysis")
        self.geometry("1540x940")
        self.minsize(1150, 720)
        self.materials = load_material_library(self.base_dir)
        self.detectors = load_detector_definitions(self.base_dir)
        self.source_presets = load_source_presets(self.base_dir)
        if not self.materials:
            messagebox.showerror("No spectra", "No absorbance spectra could be loaded.", parent=self)
            self.destroy()
            return
        self.responsivity_curve = None
        self.source_curve = None
        self.optical_terms = []
        self.result = None
        self._build_ui()
        self.apply_detector_defaults()
        self._material_changed()
        self.update_analysis()

    # ------------------------------------------------------------------ layout
    def _build_ui(self):
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
        controls = ctk.CTkScrollableFrame(self, width=360, label_text="Analysis Settings")
        controls.grid(row=0, column=0, sticky="nsew", padx=(10, 5), pady=10)
        view = ctk.CTkFrame(self)
        view.grid(row=0, column=1, sticky="nsew", padx=(5, 10), pady=10)
        self.controls = controls

        def heading(text):
            ctk.CTkLabel(controls, text=text, anchor="w", font=("Arial", 15, "bold"),
                         text_color="#7dd3fc").pack(fill="x", padx=8, pady=(16, 0))

        def option(label, values, default, command=None):
            ctk.CTkLabel(controls, text=label, anchor="w", font=("Arial", 13, "bold")).pack(fill="x", padx=8, pady=(8, 2))
            variable = ctk.StringVar(value=default)
            ctk.CTkOptionMenu(controls, variable=variable, values=list(values), command=command).pack(fill="x", padx=8)
            return variable

        def entry(label, default=""):
            ctk.CTkLabel(controls, text=label, anchor="w").pack(fill="x", padx=8, pady=(7, 1))
            variable = ctk.StringVar(value=default)
            ctk.CTkEntry(controls, textvariable=variable).pack(fill="x", padx=8)
            return variable

        def note(color="#a9c8e8"):
            label = ctk.CTkLabel(controls, text="", text_color=color, wraplength=330, justify="left", anchor="w")
            label.pack(fill="x", padx=8, pady=(3, 0))
            return label

        def checkbox(text):
            variable = ctk.BooleanVar(value=False)
            ctk.CTkCheckBox(controls, text=text, variable=variable).pack(anchor="w", padx=8, pady=(8, 2))
            return variable

        def browse_row(text, command):
            ctk.CTkButton(controls, text=text, command=command, height=28).pack(fill="x", padx=8, pady=(6, 0))
            return note()

        heading("1. Material")
        self.material_var = option("Material spectrum", self.materials, next(iter(self.materials)), self._material_changed)
        self.material_note = note()
        self.convention_var = option("Absorbance convention", CONVENTIONS, "From material data")
        self.reference_thickness_var = entry("Reference-sample thickness x_ref (µm, blank = material data)")

        heading("2. Thickness")
        self.thickness_mode_var = option("Thickness input", THICKNESS_MODES, "Layers × thickness ratio (x_ref unknown)")
        range_frame = ctk.CTkFrame(controls, fg_color="transparent")
        range_frame.pack(fill="x", padx=8, pady=(8, 0))
        self.start_var, self.end_var, self.step_var = (ctk.StringVar(value=v) for v in ("1", "4", "1"))
        for label, variable in [("Start", self.start_var), ("End", self.end_var), ("Step", self.step_var)]:
            part = ctk.CTkFrame(range_frame, fg_color="transparent")
            part.pack(side="left", expand=True, fill="x", padx=2)
            ctk.CTkLabel(part, text=label).pack()
            ctk.CTkEntry(part, textvariable=variable, width=80).pack(fill="x")
        self.layer_thickness_var = entry("Layer thickness t_layer (µm, physical mode)")
        self.ratio_var = entry("Assumed t_layer / x_ref ratio (ratio mode)", "1")
        self.inspect_var = entry("Value to inspect (must be one of the simulated values)", "3")

        heading("3. Detector")
        self.detector_var = option("Detector dataset (defaults + measured voltages)", self.detectors,
                                   next(iter(self.detectors)), lambda _v: self._detector_changed())
        ctk.CTkButton(controls, text="Apply detector defaults (V0, V_dark, band)", command=self.apply_detector_defaults,
                      height=28).pack(fill="x", padx=8, pady=(6, 0))
        self.detector_note = note()
        self.v0_var = entry("No-film voltage V0 (mV)")
        self.vdark_var = entry("Dark voltage V_dark (mV, 0 unless measured)", "0")
        self.detector_mode_var = option("Detector weighting", DETECTOR_MODES, "Flat band approximation")
        self.band_min_var = entry("Wavelength window min (nm; required for flat band)")
        self.band_max_var = entry("Wavelength window max (nm; required for flat band)")
        self.responsivity_note = browse_row("Load detector responsivity CSV…", self.load_responsivity)

        heading("4. Source")
        self.source_mode_var = option("Source weighting", SOURCE_MODES, "Flat (no source weighting)")
        self.preset_var = option("Blackbody preset", [CUSTOM_TEMPERATURE, *self.source_presets], CUSTOM_TEMPERATURE,
                                 self._preset_selected)
        self.temperature_var = entry("Blackbody temperature (K)")
        self.source_min_var = entry("Source emission min (nm, optional)")
        self.source_max_var = entry("Source emission max (nm, optional)")
        self.source_note = browse_row("Load measured source spectrum CSV…", self.load_source)
        ctk.CTkButton(controls, text="Add filter/window transmission CSV…", command=self.add_optical_term,
                      height=28).pack(fill="x", padx=8, pady=(10, 0))
        ctk.CTkButton(controls, text="Clear optical terms", command=self.clear_optical_terms,
                      height=28, fg_color="#555555").pack(fill="x", padx=8, pady=(4, 0))
        self.terms_note = note()

        heading("5. Optional corrections (off by default)")
        self.baseline_var = checkbox("Subtract baseline offset (5th percentile in band)")
        self.clamp_var = checkbox("Clamp negative absorbance to zero")
        self.interface_var = checkbox("Fresnel interface loss per layer")
        self.index_var = entry("Film refractive index n (interface correction)")
        ctk.CTkLabel(controls, text=TYPICAL_INDEX_NOTE, text_color="gray", anchor="w").pack(fill="x", padx=8)

        heading("6. Comparison")
        self.compare_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(controls, text="Compare with measured detector voltages", variable=self.compare_var).pack(anchor="w", padx=8, pady=6)
        ctk.CTkButton(controls, text="UPDATE ANALYSIS", command=self.update_analysis, fg_color="#00aa00",
                      hover_color="#008800", font=("Arial", 14, "bold")).pack(fill="x", padx=8, pady=(12, 5))
        ctk.CTkButton(controls, text="EXPORT EXCEL", command=self.export_excel, fg_color="#1f6f43",
                      hover_color="#278b55", font=("Arial", 14, "bold")).pack(fill="x", padx=8, pady=(5, 12))

        view.grid_columnconfigure(0, weight=1)
        view.grid_rowconfigure(3, weight=1)
        self.status_label = ctk.CTkLabel(view, text="", font=("Arial", 14, "bold"), anchor="w", justify="left",
                                         corner_radius=6, wraplength=1050)
        self.status_label.grid(row=0, column=0, sticky="ew", padx=10, pady=(10, 0))
        self.headline_label = ctk.CTkLabel(view, text="", font=("Arial", 17, "bold"), anchor="w", justify="left")
        self.headline_label.grid(row=1, column=0, sticky="ew", padx=12, pady=(6, 0))
        graph_bar = ctk.CTkFrame(view, fg_color="transparent")
        graph_bar.grid(row=2, column=0, sticky="ew", padx=8, pady=(6, 0))
        ctk.CTkLabel(graph_bar, text="Graph:", font=("Arial", 13, "bold")).pack(side="left", padx=(4, 6))
        self.graph_var = ctk.StringVar(value=GRAPH_COMPARE)
        ctk.CTkOptionMenu(graph_bar, variable=self.graph_var, values=GRAPH_TYPES, width=360,
                          command=lambda _v: self.draw_graph()).pack(side="left")
        self.fig, self.ax = plt.subplots(figsize=(10, 5.5), facecolor="#2b2b2b")
        self.canvas = FigureCanvasTkAgg(self.fig, master=view)
        self.canvas.get_tk_widget().grid(row=3, column=0, sticky="nsew", padx=8, pady=8)
        self.results_text = ctk.CTkTextbox(view, height=210, font=("Consolas", 12), wrap="none")
        self.results_text.grid(row=4, column=0, sticky="ew", padx=8, pady=(0, 8))

    # ------------------------------------------------------- explicit actions
    def _material_changed(self, _value=None):
        material = self.materials[self.material_var.get()]
        x_ref = material["reference_thickness_um"]
        self.material_note.configure(text=(
            f"{material['source']}\nConvention: {material['convention']}   "
            f"Range: {material['wavelength_nm'][0]:.0f}-{material['wavelength_nm'][-1]:.0f} nm\n"
            f"x_ref: {'UNKNOWN' if x_ref is None else f'{x_ref:g} µm'}\n{material['note']}"))

    def _detector_changed(self):
        detector = self.detectors[self.detector_var.get()]
        differs = []
        if str(detector.get("no_film_voltage_mv")) != self.v0_var.get().strip():
            differs.append(f"V0={detector.get('no_film_voltage_mv')} mV")
        text = detector.get("responsivity_note", "")
        if differs:
            text = (f"Current settings were NOT changed. This detector's defaults: {', '.join(differs)} — "
                    f"click 'Apply detector defaults' to use them.\n" + text)
        self.detector_note.configure(text=text)

    def apply_detector_defaults(self):
        detector = self.detectors[self.detector_var.get()]
        self.v0_var.set(str(detector.get("no_film_voltage_mv", "")))
        self.vdark_var.set(str(detector.get("dark_voltage_mv") or 0))
        self.band_min_var.set(str(detector.get("flat_band_min_nm", "")))
        self.band_max_var.set(str(detector.get("flat_band_max_nm", "")))
        if detector.get("responsivity_csv"):
            self._set_responsivity(os.path.join(self.base_dir, detector["responsivity_csv"]))
        self._detector_changed()

    def _preset_selected(self, label):
        preset = self.source_presets.get(label)
        if preset is None:
            return
        self.source_mode_var.set("Blackbody approximation")
        self.temperature_var.set(f"{preset['temperature_k']:g}")
        self.source_min_var.set("" if preset["min_nm"] is None else f"{preset['min_nm']:g}")
        self.source_max_var.set("" if preset["max_nm"] is None else f"{preset['max_nm']:g}")

    def _ask_csv(self, title):
        return filedialog.askopenfilename(parent=self, title=title,
                                          filetypes=[("CSV / text", "*.csv *.txt"), ("All files", "*.*")])

    def _set_responsivity(self, path):
        self.responsivity_curve = load_spectral_curve_csv(path, "responsivity")
        self.detector_mode_var.set("Loaded detector responsivity (CSV)")
        self.responsivity_note.configure(text=self.responsivity_curve["note"])

    def load_responsivity(self):
        path = self._ask_csv("Detector responsivity CSV (wavelength, responsivity)")
        if path:
            try:
                self._set_responsivity(path)
            except Exception as exc:
                messagebox.showerror("Responsivity load failed", str(exc), parent=self)

    def load_source(self):
        path = self._ask_csv("Source spectrum CSV (wavelength, relative power)")
        if not path:
            return
        try:
            self.source_curve = load_spectral_curve_csv(path, "source")
            self.source_mode_var.set("Measured source spectrum (CSV)")
            self.source_note.configure(text=self.source_curve["note"])
        except Exception as exc:
            messagebox.showerror("Source load failed", str(exc), parent=self)

    def add_optical_term(self):
        path = self._ask_csv("Filter/window transmission CSV (wavelength, fraction 0-1)")
        if not path:
            return
        try:
            self.optical_terms.append(load_spectral_curve_csv(path, "transmission"))
            self.terms_note.configure(text="\n".join(term["note"] for term in self.optical_terms))
        except Exception as exc:
            messagebox.showerror("Optical term load failed", str(exc), parent=self)

    def clear_optical_terms(self):
        self.optical_terms = []
        self.terms_note.configure(text="")

    # ------------------------------------------------------------- analysis
    def _values(self):
        start = _required_float(self.start_var.get(), "Start")
        end = _required_float(self.end_var.get(), "End")
        step = _required_float(self.step_var.get(), "Step")
        if step <= 0 or end < start:
            raise ValueError("Thickness range requires end ≥ start and step > 0.")
        count = int(math.floor((end - start) / step + 1e-9)) + 1
        if count > 1000:
            raise ValueError("Thickness range is limited to 1000 points.")
        return [round(start + i * step, 10) for i in range(count)]

    def build_config(self):
        material_name = self.material_var.get()
        detector = self.detectors[self.detector_var.get()]
        measured = None
        if self.compare_var.get():
            measured = {int(k): float(v) for k, v in detector.get("materials", {}).get(material_name, {}).items()}
        return BroadbandConfig(
            material=material_name,
            values=self._values(),
            no_film_voltage_mv=_required_float(self.v0_var.get(), "V0"),
            dark_voltage_mv=_optional_float(self.vdark_var.get(), "V_dark") or 0.0,
            thickness_mode=THICKNESS_MODES[self.thickness_mode_var.get()],
            reference_thickness_um=_optional_float(self.reference_thickness_var.get(), "x_ref"),
            layer_thickness_um=_optional_float(self.layer_thickness_var.get(), "Layer thickness"),
            layer_to_reference_ratio=_optional_float(self.ratio_var.get(), "Thickness ratio"),
            wavelength_min_nm=_optional_float(self.band_min_var.get(), "Window min"),
            wavelength_max_nm=_optional_float(self.band_max_var.get(), "Window max"),
            detector_name=detector.get("display_name", self.detector_var.get()),
            detector_mode=DETECTOR_MODES[self.detector_mode_var.get()],
            detector_responsivity=self.responsivity_curve,
            source_mode=SOURCE_MODES[self.source_mode_var.get()],
            source_label="" if self.preset_var.get() == CUSTOM_TEMPERATURE else self.preset_var.get(),
            source_temperature_k=_optional_float(self.temperature_var.get(), "Blackbody temperature"),
            source_min_nm=_optional_float(self.source_min_var.get(), "Source min"),
            source_max_nm=_optional_float(self.source_max_var.get(), "Source max"),
            source_spectrum=self.source_curve,
            optical_terms=list(self.optical_terms),
            absorbance_convention=CONVENTIONS[self.convention_var.get()],
            baseline_correction=self.baseline_var.get(),
            clamp_negative=self.clamp_var.get(),
            interface_correction=self.interface_var.get(),
            refractive_index=_optional_float(self.index_var.get(), "Refractive index"),
            measured_voltages_by_layer=measured,
            inspect_value=_optional_float(self.inspect_var.get(), "Inspected value"),
        )

    def update_analysis(self):
        try:
            self.result = run_broadband_analysis(self.build_config(), self.materials)
        except Exception as exc:
            messagebox.showerror("Analysis error", str(exc), parent=self)
            return
        self._update_status()
        self._update_results()
        self.draw_graph()

    def _thickness_label(self, row):
        if row["thickness_um"] is not None:
            return f"{row['layer_count']} layer(s) = {row['thickness_um']:g} µm"
        if row["layer_count"] is not None:
            return f"{row['layer_count']} layer(s) (x/x_ref = {row['thickness_scale']:g})"
        return f"x/x_ref = {row['thickness_scale']:g}"

    def _x_axis(self):
        mode = self.result["thickness"]["mode"]
        rows = self.result["rows"]
        if mode == LAYERS_PHYSICAL:
            return [r["thickness_um"] for r in rows], "Thickness x (µm)"
        if mode == LAYERS_RATIO:
            return [r["layer_count"] for r in rows], "Layer count (assumed t_layer/x_ref ratio)"
        return [r["thickness_scale"] for r in rows], "Equivalent reference samples x/x_ref"

    def _update_status(self):
        result = self.result
        coverage = result["coverage"]
        corrections = result["corrections"]
        lines = [f"Spectral coverage: {coverage['fraction']:.1%} of detector/source weight "
                 f"(weight {coverage['weight_range_nm'][0]:.0f}-{coverage['weight_range_nm'][1]:.0f} nm, "
                 f"material {coverage['covered_range_nm'][0]:.0f}-{coverage['covered_range_nm'][1]:.0f} nm)"]
        lines += result["warnings"]
        applied = []
        if corrections["baseline_enabled"]:
            applied.append(f"baseline offset {corrections['baseline_offset']:+.4f} A subtracted")
        if corrections["clamp_enabled"]:
            applied.append(f"{corrections['clamped_points']} points ({corrections['clamped_percent']:.1f}%) clamped to 0")
        if corrections["interface_enabled"]:
            applied.append(f"interface T = {corrections['interface_per_layer']:.4f} per layer (n = {corrections['refractive_index']:g}); "
                           + corrections["interface_assumptions"])
        lines.append("Corrections applied: " + ("; ".join(applied) if applied else "none (raw forward model)"))
        color = "#7f1d1d" if coverage["partial"] else ("#78350f" if result["warnings"] else "#14532d")
        self.status_label.configure(text="\n".join(lines), fg_color=color)

        index = result["inspected_index"]
        if index is None:
            self.headline_label.configure(text="")
            return
        row = result["rows"][index]
        tag = "  [PARTIAL COVERAGE — not fully valid]" if row["partial_coverage"] else ""
        text = f"{result['material']['name']}, {self._thickness_label(row)}:   Raw prediction {row['raw_predicted_voltage_mv']:.1f} mV"
        if corrections["any_enabled"]:
            text += f"   |   Corrected {row['predicted_voltage_mv']:.1f} mV"
        if row["measured_voltage_mv"] is not None:
            text += f"   |   Measured {row['measured_voltage_mv']:.1f} mV   |   Error raw {row['raw_error_mv']:+.1f} mV / {row['raw_percent_error']:+.1f}%"
            if corrections["any_enabled"]:
                text += f", corrected {row['error_mv']:+.1f} mV / {row['percent_error']:+.1f}%"
        self.headline_label.configure(text=text + tag)

    def _update_results(self):
        result = self.result
        fmt = lambda v, spec: "—" if v is None else format(v, spec)
        lines = [f"Thickness: {result['thickness']['description']}",
                 f"Detector:  {result['weighting']['detector_description']}",
                 f"Source:    {result['weighting']['source_description']}",
                 f"V0 = {result['voltages']['no_film_voltage_mv']:g} mV, V_dark = {result['voltages']['dark_voltage_mv']:g} mV;  "
                 f"Vpred = V_dark + (V0 - V_dark)·T_eff",
                 "",
                 f"{'Layers':>6} {'x (µm)':>8} {'x/x_ref':>8} {'Meas mV':>9} {'Raw mV':>9} {'Corr mV':>9} "
                 f"{'T_spec':>7} {'T_int':>7} {'T_eff':>7} {'Err mV':>8} {'Err %':>7} {'Cover':>6}"]
        for row in result["rows"]:
            lines.append(f"{fmt(row['layer_count'], 'd'):>6} {fmt(row['thickness_um'], '.1f'):>8} {row['thickness_scale']:>8.3f} "
                         f"{fmt(row['measured_voltage_mv'], '.1f'):>9} {row['raw_predicted_voltage_mv']:>9.1f} "
                         f"{row['predicted_voltage_mv']:>9.1f} {row['spectral_transmission']:>7.4f} "
                         f"{row['interface_transmission']:>7.4f} {row['predicted_transmission']:>7.4f} "
                         f"{fmt(row['error_mv'], '+.1f'):>8} {fmt(row['percent_error'], '+.1f'):>7} "
                         f"{row['coverage_fraction']:>6.1%}")
        metric_sets = [("Raw", result["raw_metrics"])]
        if result["corrections"]["any_enabled"]:
            metric_sets.append(("Corrected", result["metrics"]))
        for label, metrics in metric_sets:
            if metrics["n"]:
                lines.append(f"{label:>9} agreement (n={metrics['n']}): MAE={metrics['mae_mv']:.2f} mV  "
                             f"RMSE={metrics['rmse_mv']:.2f} mV  MAPE={fmt(metrics['mape_percent'], '.1f')}%  "
                             f"R²={fmt(metrics['r_squared'], '.4f')}"
                             + ("  [PARTIAL COVERAGE]" if result["coverage"]["partial"] else ""))
        self.results_text.configure(state="normal")
        self.results_text.delete("1.0", "end")
        self.results_text.insert("1.0", "\n".join(lines))
        self.results_text.configure(state="disabled")

    # ------------------------------------------------------------- plotting
    def _style_axis(self, title, xlabel, ylabel):
        self.ax.clear()
        self.ax.set_autoscale_on(True)
        self.ax.set_facecolor("#2b2b2b")
        self.ax.tick_params(colors="white")
        for spine in self.ax.spines.values():
            spine.set_color("gray")
        suffix = "  [PARTIAL COVERAGE]" if self.result["coverage"]["partial"] else ""
        self.ax.set_title(title + suffix, color="#fca5a5" if suffix else "white")
        self.ax.set_xlabel(xlabel, color="white")
        self.ax.set_ylabel(ylabel, color="white")
        self.ax.grid(True, alpha=0.2)

    def _legend(self, title=None):
        if self.ax.get_legend_handles_labels()[0]:
            self.ax.legend(title=title, facecolor="#333333", edgecolor="gray", labelcolor="white")

    def _shade_coverage(self):
        low, high = self.result["coverage"]["covered_range_nm"]
        self.ax.axvspan(low, high, color="#64748b", alpha=0.12, label="Material spectrum coverage")

    def draw_graph(self):
        if self.result is None:
            return
        result, graph = self.result, self.graph_var.get()
        grid, rows = result["grid"], result["rows"]
        wl = grid["wavelength_nm"]
        xs, xlabel = self._x_axis()
        compared = [(x, r) for x, r in zip(xs, rows) if r["measured_voltage_mv"] is not None]
        index = result["inspected_index"] if result["inspected_index"] is not None else 0
        corrected = result["corrections"]["any_enabled"]

        if graph == GRAPH_ABSORBANCE:
            self._style_axis(f"{result['material']['name']} reference absorbance ({result['material']['convention']})",
                             "Wavelength (nm)", "Absorbance")
            self.ax.plot(wl, grid["raw_absorbance"], color="#94a3b8" if corrected else "#38bdf8", linewidth=1, label="Raw")
            if corrected:
                self.ax.plot(wl, grid["absorbance_used"], color="#38bdf8", label="After enabled corrections")
            self.ax.axhline(0, color="gray", linewidth=0.8)
            self._legend()
        elif graph == GRAPH_TRANSMISSION:
            self._style_axis("Spectral transmission T(λ, x)  (interface factor not included)", "Wavelength (nm)", "T")
            show = range(len(rows)) if len(rows) <= len(SERIES_COLORS) else [index]
            for color_index, i in enumerate(show):
                self.ax.plot(wl, result["spectra"][i]["transmission"], color=SERIES_COLORS[color_index % len(SERIES_COLORS)],
                             linewidth=2.2 if i == index else 1.1, label=self._thickness_label(rows[i]))
            self._legend()
        elif graph == GRAPH_WEIGHTING:
            self._style_axis("Detector / source weighting", "Wavelength (nm)", "Normalized value")
            self._shade_coverage()
            for values, label, color in [(grid["source"], "Source S(λ)", "#fbbf24"),
                                         (grid["responsivity"], "Detector R(λ)", "#38bdf8")]:
                if values is not None:
                    self.ax.plot(wl, values / (np.max(values) or 1), label=label, color=color)
            for i, values in enumerate(grid["optical_terms"]):
                self.ax.plot(wl, values, label=f"Optical term {i + 1}", color="#a78bfa")
            self.ax.plot(wl, grid["weight"] / (np.max(grid["weight"]) or 1), label="W(λ) normalized", color="#f472b6", linewidth=2)
            self.ax.set_ylim(0, 1.08)
            self._legend()
        elif graph == GRAPH_CONTRIBUTION:
            spectral = result["spectra"][index]
            peak = np.max(grid["weight"]) or 1.0
            self._style_axis(f"Weighted contribution at {self._thickness_label(rows[index])}", "Wavelength (nm)",
                             "Normalized contribution")
            self.ax.fill_between(wl, grid["weight"] / peak, color="#38bdf8", alpha=0.2, label="No film: W(λ)")
            self.ax.fill_between(wl, np.nan_to_num(spectral["weighted_transmitted"]) / peak, color="#fb7185", alpha=0.6,
                                 label="With film: W(λ)·T(λ)")
            self.ax.text(0.01, 0.03, f"T_spectral = {spectral['spectral_transmission']:.1%}   coverage = "
                         f"{spectral['coverage_fraction']:.1%}", transform=self.ax.transAxes, color="white")
            self._legend()
        elif graph in (GRAPH_VOLTAGE, GRAPH_COMPARE):
            self._style_axis(graph, xlabel, "Detector voltage (mV)")
            self.ax.plot(xs, [r["raw_predicted_voltage_mv"] for r in rows], marker="o", color="#34d399", label="Raw prediction")
            if corrected:
                self.ax.plot(xs, [r["predicted_voltage_mv"] for r in rows], marker="s", color="#38bdf8", label="Corrected prediction")
            if result["coverage"]["partial"]:
                self.ax.fill_between(xs, [r["predicted_voltage_bounds_mv"][0] for r in rows],
                                     [r["predicted_voltage_bounds_mv"][1] for r in rows], color="#f87171", alpha=0.15,
                                     label="Bounds (uncovered weight T=0…1)")
            if graph == GRAPH_COMPARE and compared:
                self.ax.scatter([x for x, _ in compared], [r["measured_voltage_mv"] for _, r in compared],
                                color="#fb7185", s=60, zorder=3, label="Measured")
            self._legend()
        elif graph == GRAPH_PARITY:
            self._style_axis(graph, "Measured voltage (mV)", "Predicted voltage (mV)")
            if compared:
                measured = [r["measured_voltage_mv"] for _, r in compared]
                keys = [("raw_predicted_voltage_mv", "Raw", "#34d399")]
                if corrected:
                    keys.append(("predicted_voltage_mv", "Corrected", "#38bdf8"))
                values = measured + [r[key] for key, _, _ in keys for _, r in compared]
                low, high = min(values) * 0.95, max(values) * 1.05
                self.ax.plot([low, high], [low, high], color="gray", linestyle="--", label="1:1")
                for key, label, color in keys:
                    self.ax.scatter(measured, [r[key] for _, r in compared], color=color, s=60, zorder=3, label=label)
                self._legend()
        elif graph == GRAPH_ERROR:
            self._style_axis(graph, xlabel, "Error (%)  [+ = over-predicts]")
            if compared:
                width = 0.35 * (min(np.diff(sorted(set(xs)))) if len(set(xs)) > 1 else 1)
                self.ax.bar([x - (width / 2 if corrected else 0) for x, _ in compared],
                            [r["raw_percent_error"] for _, r in compared], width=width, color="#34d399", label="Raw")
                if corrected:
                    self.ax.bar([x + width / 2 for x, _ in compared], [r["percent_error"] for _, r in compared],
                                width=width, color="#38bdf8", label="Corrected")
                self.ax.axhline(0, color="gray")
                self._legend()
        elif graph == GRAPH_EFFECTIVE_T:
            self._style_axis(graph, xlabel, "T_eff = (V - V_dark)/(V0 - V_dark)")
            self.ax.plot(xs, [r["raw_effective_transmission"] for r in rows], marker="o", color="#34d399", label="Raw prediction")
            if corrected:
                self.ax.plot(xs, [r["predicted_transmission"] for r in rows], marker="s", color="#38bdf8", label="Corrected")
            if compared:
                self.ax.scatter([x for x, _ in compared], [r["measured_transmission"] for _, r in compared],
                                color="#fb7185", s=60, zorder=3, label="Measured")
            self._legend()
        else:
            self._style_axis(graph, xlabel, "ln(1/T_eff)")
            self.ax.plot(xs, [np.nan if r["predicted_attenuation_natural"] is None else r["predicted_attenuation_natural"]
                              for r in rows], marker="o", color="#a78bfa", label="Predicted (corrections as enabled)")
            if compared:
                self.ax.scatter([x for x, _ in compared],
                                [np.nan if r["measured_attenuation_natural"] is None else r["measured_attenuation_natural"]
                                 for _, r in compared], color="#fb7185", s=60, zorder=3, label="Measured")
            self._legend()
        self.fig.tight_layout()
        self.canvas.draw_idle()

    # --------------------------------------------------------------- export
    def experimental_rows(self):
        rows = []
        for detector in self.detectors.values():
            for material, measurements in detector.get("materials", {}).items():
                for layer, voltage in measurements.items():
                    rows.append({"material": material, "detector": detector["display_name"], "layer_count": int(layer),
                                 "measured_voltage_mv": voltage, "no_film_voltage_mv": detector["no_film_voltage_mv"],
                                 "dark_voltage_mv": detector.get("dark_voltage_mv")})
        return rows

    def export_excel(self):
        if self.result is None:
            return
        default_name = f"broadband_{self.material_var.get()}_{self.detector_var.get()}.xlsx"
        default_name = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in default_name)
        path = filedialog.asksaveasfilename(parent=self, title="Export broadband analysis", initialfile=default_name,
                                            defaultextension=".xlsx", filetypes=[("Excel workbook", "*.xlsx")])
        if not path:
            return
        try:
            saved = export_broadband_workbook(path, self.result, self.experimental_rows())
            messagebox.showinfo("Export complete", f"Workbook saved to:\n{saved}", parent=self)
        except Exception as exc:
            messagebox.showerror("Export failed", str(exc), parent=self)
