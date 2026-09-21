"""CustomTkinter window for broadband detector analysis and export.

The visible panel holds only the five things a normal run needs: material,
detector, mode, one value and the buttons.  Everything else lives behind
"Advanced settings", with defaults that make the basic panel sufficient.

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
from BroadbandPipeline import (BASELINE_AUTO, BASELINE_OFF, BASELINE_ON, DETECTOR_FLAT_BAND,
                               DETECTOR_RESPONSIVITY, LAYERS_PHYSICAL, LAYERS_RATIO, REFERENCE_MULTIPLIER,
                               SOURCE_BLACKBODY, SOURCE_FLAT, SOURCE_MEASURED, BroadbandConfig,
                               run_broadband_analysis)
from BroadbandSheetUI import choose_measured_points
from SpectralData import (load_detector_definitions, load_material_library, load_source_presets,
                          load_spectral_curve_csv)

# Measured points are drawn on every thickness graph, so there is one voltage
# graph rather than a "predicted" and a "measured vs predicted" version of it.
GRAPH_VOLTAGE = "Voltage vs thickness"
GRAPH_EFFECTIVE_T = "Transmission vs thickness"
GRAPH_ERROR = "Prediction error (%)"
GRAPH_PARITY = "Parity: measured vs predicted"
GRAPH_ATTENUATION = "Attenuation ln(1/T) vs thickness"
GRAPH_ABSORBANCE = "Material absorbance spectrum"
GRAPH_TRANSMISSION = "Transmission spectrum by thickness"
GRAPH_WEIGHTING = "Detector + source weighting"
GRAPH_CONTRIBUTION = "Where the signal comes from"
GRAPH_TYPES = [GRAPH_VOLTAGE, GRAPH_EFFECTIVE_T, GRAPH_ERROR, GRAPH_PARITY, GRAPH_ATTENUATION,
               GRAPH_ABSORBANCE, GRAPH_TRANSMISSION, GRAPH_WEIGHTING, GRAPH_CONTRIBUTION]

MEASURED_BUILTIN = "Built-in detector data"
MEASURED_FILE = "From an Excel/CSV file\u2026"
MEASURED_SOURCES = [MEASURED_BUILTIN, MEASURED_FILE]

MODE_PREDICT = "Predict voltage"
MODE_ESTIMATE = "Estimate thickness"
MODES = [MODE_PREDICT, MODE_ESTIMATE]
MODE_HINTS = {
    MODE_PREDICT: "Enter a film thickness; the model predicts the detector voltage.",
    MODE_ESTIMATE: "Enter a measured detector voltage; the same model is inverted for the thickness.",
}

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
BASELINE_MODES = {"Automatic (only if the data needs it)": BASELINE_AUTO,
                  "Always subtract the baseline": BASELINE_ON,
                  "Never correct the baseline": BASELINE_OFF}
CONVENTIONS = {"From material data": None, "Override: base-10 absorbance": "base10",
               "Override: natural-log attenuation": "natural"}
CUSTOM_TEMPERATURE = "Custom temperature"
TYPICAL_INDEX_NOTE = "Typical mid-IR n: PE ≈ 1.51, EVOH ≈ 1.52, Nylon ≈ 1.53"
SERIES_COLORS = ["#38bdf8", "#34d399", "#fbbf24", "#fb7185", "#a78bfa", "#f472b6", "#2dd4bf", "#fb923c"]

# The graph always shows a sweep around the value being analysed.
MIN_SWEEP_POINTS = 4
MAX_SWEEP_POINTS = 20


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
        self.geometry("1440x900")
        self.minsize(1080, 700)
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
        self.imported_measurements = None      # {"points", "source", "no_film_voltage_mv", ...}
        self.result = None
        self.advanced_visible = False
        self.details_visible = False
        self._build_ui()
        self.apply_detector_defaults()
        self._material_changed()
        self.update_analysis()

    # ------------------------------------------------------------------ layout
    def _build_ui(self):
        self.grid_columnconfigure(2, weight=1)
        self.grid_rowconfigure(0, weight=1)
        basics = ctk.CTkFrame(self, width=320)
        basics.grid(row=0, column=0, sticky="nsew", padx=(10, 5), pady=10)
        basics.grid_propagate(False)
        self.advanced = ctk.CTkScrollableFrame(self, width=330, label_text="Advanced settings")
        view = ctk.CTkFrame(self)
        view.grid(row=0, column=2, sticky="nsew", padx=(5, 10), pady=10)

        def option(parent, label, values, default, command=None):
            ctk.CTkLabel(parent, text=label, anchor="w", font=("Arial", 13, "bold")).pack(fill="x", padx=10, pady=(10, 2))
            variable = ctk.StringVar(value=default)
            ctk.CTkOptionMenu(parent, variable=variable, values=list(values), command=command).pack(fill="x", padx=10)
            return variable

        def entry(parent, label, default=""):
            ctk.CTkLabel(parent, text=label, anchor="w").pack(fill="x", padx=10, pady=(7, 1))
            variable = ctk.StringVar(value=default)
            ctk.CTkEntry(parent, textvariable=variable).pack(fill="x", padx=10)
            return variable

        def note(parent, color="#a9c8e8", size=11):
            label = ctk.CTkLabel(parent, text="", text_color=color, font=("Arial", size), wraplength=296,
                                 justify="left", anchor="w")
            label.pack(fill="x", padx=10, pady=(4, 0))
            return label

        def heading(parent, text):
            ctk.CTkLabel(parent, text=text, anchor="w", font=("Arial", 14, "bold"),
                         text_color="#7dd3fc").pack(fill="x", padx=10, pady=(14, 0))

        def checkbox(parent, text, default=False):
            variable = ctk.BooleanVar(value=default)
            ctk.CTkCheckBox(parent, text=text, variable=variable).pack(anchor="w", padx=10, pady=(8, 2))
            return variable

        def browse_row(parent, text, command):
            ctk.CTkButton(parent, text=text, command=command, height=28).pack(fill="x", padx=10, pady=(6, 0))
            return note(parent)

        # ------------------------------------------------------------ basics
        ctk.CTkLabel(basics, text="Broadband analysis", anchor="w",
                     font=("Arial", 17, "bold")).pack(fill="x", padx=10, pady=(12, 0))
        self.material_var = option(basics, "1. Material", self.materials, next(iter(self.materials)),
                                   self._material_changed)
        self.detector_var = option(basics, "2. Detector", self._detector_labels(), self._detector_labels()[0],
                                   lambda _v: self._detector_changed())
        self.measured_source_var = option(basics, "Measured data", MEASURED_SOURCES, MEASURED_BUILTIN,
                                          self._measured_source_changed)
        self.measured_note = note(basics, "#94a3b8")
        ctk.CTkLabel(basics, text="3. Mode", anchor="w", font=("Arial", 13, "bold")).pack(fill="x", padx=10, pady=(12, 2))
        self.mode_var = ctk.StringVar(value=MODE_PREDICT)
        ctk.CTkSegmentedButton(basics, values=MODES, variable=self.mode_var,
                               command=lambda _v: self._mode_changed()).pack(fill="x", padx=10)
        self.mode_hint = note(basics, "#94a3b8")
        self.value_label = ctk.CTkLabel(basics, text="", anchor="w", font=("Arial", 13, "bold"))
        self.value_label.pack(fill="x", padx=10, pady=(12, 2))
        self.value_var = ctk.StringVar(value="3")
        value_entry = ctk.CTkEntry(basics, textvariable=self.value_var, font=("Arial", 15))
        value_entry.pack(fill="x", padx=10)
        value_entry.bind("<Return>", lambda _e: self.update_analysis())
        ctk.CTkButton(basics, text="ANALYZE", command=self.update_analysis, fg_color="#00aa00",
                      hover_color="#008800", height=40, font=("Arial", 15, "bold")).pack(fill="x", padx=10, pady=(16, 6))
        ctk.CTkButton(basics, text="Export to Excel…", command=self.export_excel, fg_color="#1f6f43",
                      hover_color="#278b55", height=34, font=("Arial", 13, "bold")).pack(fill="x", padx=10)
        self.advanced_button = ctk.CTkButton(basics, text="Advanced settings  ▸", command=self.toggle_advanced,
                                             fg_color="#3f3f46", hover_color="#52525b", height=30)
        self.advanced_button.pack(fill="x", padx=10, pady=(16, 4))
        self.setup_note = note(basics, "#94a3b8")
        self.material_note = note(basics, "#6b7280", size=10)

        # ---------------------------------------------------------- advanced
        heading(self.advanced, "Material")
        self.convention_var = option(self.advanced, "Absorbance convention", CONVENTIONS, "From material data")
        self.reference_thickness_var = entry(self.advanced, "Reference thickness x_ref (µm, blank = material data)")

        heading(self.advanced, "Thickness definition")
        self.thickness_mode_var = option(self.advanced, "Thickness input", THICKNESS_MODES,
                                         "Layers × thickness ratio (x_ref unknown)", lambda _v: self._mode_changed())
        self.layer_thickness_var = entry(self.advanced, "Layer thickness t_layer (µm, physical mode)")
        self.ratio_var = entry(self.advanced, "Assumed t_layer / x_ref ratio (ratio mode)", "1")

        heading(self.advanced, "Detector voltages & band")
        self.detector_note = note(self.advanced)
        self.v0_var = entry(self.advanced, "No-film voltage V0 (mV)")
        self.vdark_var = entry(self.advanced, "Dark voltage V_dark (mV, 0 unless measured)", "0")
        ctk.CTkButton(self.advanced, text="Reset to detector defaults", command=self.apply_detector_defaults,
                      height=28).pack(fill="x", padx=10, pady=(6, 0))
        self.detector_mode_var = option(self.advanced, "Detector weighting", DETECTOR_MODES, "Flat band approximation")
        self.band_min_var = entry(self.advanced, "Wavelength window min (nm)")
        self.band_max_var = entry(self.advanced, "Wavelength window max (nm)")
        self.responsivity_note = browse_row(self.advanced, "Load detector responsivity CSV…", self.load_responsivity)

        heading(self.advanced, "Source & optics")
        self.source_mode_var = option(self.advanced, "Source weighting", SOURCE_MODES, "Flat (no source weighting)")
        self.preset_var = option(self.advanced, "Blackbody preset", [CUSTOM_TEMPERATURE, *self.source_presets],
                                 CUSTOM_TEMPERATURE, self._preset_selected)
        self.temperature_var = entry(self.advanced, "Blackbody temperature (K)")
        self.source_min_var = entry(self.advanced, "Source emission min (nm, optional)")
        self.source_max_var = entry(self.advanced, "Source emission max (nm, optional)")
        self.source_note = browse_row(self.advanced, "Load measured source spectrum CSV…", self.load_source)
        ctk.CTkButton(self.advanced, text="Add filter/window transmission CSV…", command=self.add_optical_term,
                      height=28).pack(fill="x", padx=10, pady=(10, 0))
        ctk.CTkButton(self.advanced, text="Clear optical terms", command=self.clear_optical_terms, height=28,
                      fg_color="#555555").pack(fill="x", padx=10, pady=(4, 0))
        self.terms_note = note(self.advanced)

        heading(self.advanced, "Corrections")
        self.baseline_var = option(self.advanced, "Baseline offset", BASELINE_MODES,
                                   "Automatic (only if the data needs it)")
        self.baseline_note = note(self.advanced, "#94a3b8")
        self.clamp_var = checkbox(self.advanced, "Clamp negative absorbance to zero")
        self.interface_var = checkbox(self.advanced, "Fresnel interface loss per layer")
        self.index_var = entry(self.advanced, "Film refractive index n (interface correction)")
        ctk.CTkLabel(self.advanced, text=TYPICAL_INDEX_NOTE, text_color="gray", font=("Arial", 10),
                     anchor="w").pack(fill="x", padx=10)
        self.compare_var = checkbox(self.advanced, "Compare with measured detector voltages", default=True)

        # -------------------------------------------------------------- view
        view.grid_columnconfigure(0, weight=1)
        view.grid_rowconfigure(5, weight=1)
        self.headline_label = ctk.CTkLabel(view, text="", font=("Arial", 26, "bold"), anchor="w", justify="left")
        self.headline_label.grid(row=0, column=0, sticky="ew", padx=12, pady=(14, 0))
        self.subline_label = ctk.CTkLabel(view, text="", font=("Arial", 14), text_color="#cbd5e1", anchor="w",
                                          justify="left", wraplength=1000)
        self.subline_label.grid(row=1, column=0, sticky="ew", padx=12, pady=(2, 0))
        self.status_label = ctk.CTkLabel(view, text="", font=("Arial", 12), anchor="w", justify="left",
                                         corner_radius=6, wraplength=1000)
        self.status_label.grid(row=2, column=0, sticky="ew", padx=10, pady=(8, 0))
        self.explain_label = ctk.CTkLabel(view, text="", font=("Consolas", 12), anchor="w", justify="left",
                                          corner_radius=6, fg_color="#1f2937", text_color="#e2e8f0", wraplength=1000)
        graph_bar = ctk.CTkFrame(view, fg_color="transparent")
        graph_bar.grid(row=4, column=0, sticky="ew", padx=8, pady=(8, 0))
        ctk.CTkLabel(graph_bar, text="Graph:", font=("Arial", 13, "bold")).pack(side="left", padx=(4, 6))
        self.graph_var = ctk.StringVar(value=GRAPH_VOLTAGE)
        ctk.CTkOptionMenu(graph_bar, variable=self.graph_var, values=GRAPH_TYPES, width=280,
                          command=lambda _v: self.draw_graph()).pack(side="left")
        self.explain_button = ctk.CTkButton(graph_bar, text="How is this calculated?  ▾", width=200, height=28,
                                            fg_color="#3f3f46", hover_color="#52525b", command=self.toggle_explanation)
        self.explain_button.pack(side="left", padx=12)
        self.details_button = ctk.CTkButton(graph_bar, text="Numbers  ▾", width=110, height=28, fg_color="#3f3f46",
                                            hover_color="#52525b", command=self.toggle_details)
        self.details_button.pack(side="right", padx=4)
        self.fig, self.ax = plt.subplots(figsize=(10, 5.5), facecolor="#2b2b2b")
        self.canvas = FigureCanvasTkAgg(self.fig, master=view)
        self.canvas.get_tk_widget().grid(row=5, column=0, sticky="nsew", padx=8, pady=8)
        self.results_text = ctk.CTkTextbox(view, height=210, font=("Consolas", 12), wrap="none")
        self.explanation_visible = False
        self.view = view
        self._mode_changed()

    def _detector_labels(self):
        return [detector.get("display_name", key) for key, detector in self.detectors.items()]

    def _detector_key(self):
        label = self.detector_var.get()
        for key, detector in self.detectors.items():
            if detector.get("display_name", key) == label:
                return key
        return label

    # -------------------------------------------------------- show/hide panes
    def toggle_advanced(self):
        self.advanced_visible = not self.advanced_visible
        if self.advanced_visible:
            self.advanced.grid(row=0, column=1, sticky="nsew", padx=(0, 5), pady=10)
            self.advanced_button.configure(text="Advanced settings  ▾")
        else:
            self.advanced.grid_forget()
            self.advanced_button.configure(text="Advanced settings  ▸")

    def toggle_details(self):
        self.details_visible = not self.details_visible
        if self.details_visible:
            self.results_text.grid(row=6, column=0, sticky="ew", padx=8, pady=(0, 8))
            self.details_button.configure(text="Numbers  ▴")
        else:
            self.results_text.grid_forget()
            self.details_button.configure(text="Numbers  ▾")

    def toggle_explanation(self):
        self.explanation_visible = not self.explanation_visible
        if self.explanation_visible:
            self.explain_label.grid(row=3, column=0, sticky="ew", padx=10, pady=(8, 0))
            self.explain_button.configure(text="How is this calculated?  ▴")
        else:
            self.explain_label.grid_forget()
            self.explain_button.configure(text="How is this calculated?  ▾")

    # ------------------------------------------------------- explicit actions
    def _mode_changed(self, _value=None):
        multiplier = THICKNESS_MODES[self.thickness_mode_var.get()] == REFERENCE_MULTIPLIER
        if self.mode_var.get() == MODE_ESTIMATE:
            self.value_label.configure(text="4. Measured detector voltage (mV)")
        else:
            self.value_label.configure(text="4. Equivalent reference samples x/x_ref" if multiplier
                                       else "4. Number of film layers")
        self.mode_hint.configure(text=MODE_HINTS[self.mode_var.get()])
        self.value_var.set(self._default_value())

    def _default_value(self):
        """A starting value that suits the current material, detector and mode."""
        if self.mode_var.get() == MODE_ESTIMATE:
            if self._using_imported():
                return f"{self.imported_measurements['points'][0][1]:g}"
            measured = self._measured_voltages()
            if measured:
                return f"{measured[min(measured)]:g}"
            v0 = _optional_float(self.v0_var.get(), "V0")
            return f"{v0 / 2:g}" if v0 else "100"
        return "1" if THICKNESS_MODES[self.thickness_mode_var.get()] == REFERENCE_MULTIPLIER else "3"

    def _measured_voltages(self):
        """Built-in measurements for this material and detector, keyed by layer count."""
        detector = self.detectors[self._detector_key()]
        return {int(k): float(v) for k, v in detector.get("materials", {}).get(self.material_var.get(), {}).items()}

    def _using_imported(self):
        return self.measured_source_var.get() == MEASURED_FILE and self.imported_measurements is not None

    def _measured_source_changed(self, choice):
        if choice != MEASURED_FILE:
            self.imported_measurements = None
            self._update_measured_note()
            self.update_analysis()
            return
        path = filedialog.askopenfilename(
            parent=self, title="Measured data (Excel or CSV)",
            filetypes=[("Spreadsheets", "*.xlsx *.xlsm *.csv *.txt"), ("All files", "*.*")])
        if path:
            try:
                caption = "x/x_ref" if THICKNESS_MODES[self.thickness_mode_var.get()] == REFERENCE_MULTIPLIER \
                    else "film layers"
                chosen = choose_measured_points(self, path, caption)
            except Exception as exc:
                messagebox.showerror("Could not open that file", str(exc), parent=self)
                chosen = None
            if chosen:
                self.imported_measurements = chosen
                if chosen.get("no_film_voltage_mv") is not None:
                    self.v0_var.set(f"{chosen['no_film_voltage_mv']:g}")
                self._update_measured_note()
                self.update_analysis()
                return
        # Nothing usable was chosen, so stay on whatever was in use before.
        if self.imported_measurements is None:
            self.measured_source_var.set(MEASURED_BUILTIN)
        self._update_measured_note()

    def _update_measured_note(self):
        if self._using_imported():
            imported = self.imported_measurements
            text = imported["summary"]
            if imported.get("no_film_voltage_mv") is not None:
                text += f"\nNo-film voltage V0 set to {imported['no_film_voltage_mv']:g} mV from the file."
        else:
            count = len(self._measured_voltages())
            text = (f"{count} layer voltage(s) recorded for this material and detector"
                    if count else "No built-in measurements for this material and detector.")
        self.measured_note.configure(text=text)

    def _material_changed(self, _value=None):
        material = self.materials[self.material_var.get()]
        x_ref = material["reference_thickness_um"]
        self.material_note.configure(text=(
            f"Spectrum {material['wavelength_nm'][0]:.0f}-{material['wavelength_nm'][-1]:.0f} nm · "
            f"reference film thickness {'unknown' if x_ref is None else f'{x_ref:g} µm'}"))
        if self.mode_var.get() == MODE_ESTIMATE:
            self.value_var.set(self._default_value())
        self._update_measured_note()
        self._update_setup_note()

    def _detector_changed(self):
        self.apply_detector_defaults()
        if self.mode_var.get() == MODE_ESTIMATE:
            self.value_var.set(self._default_value())

    def apply_detector_defaults(self):
        """Load this detector's V0, V_dark, band and responsivity into the advanced fields."""
        detector = self.detectors[self._detector_key()]
        self.v0_var.set(str(detector.get("no_film_voltage_mv", "")))
        self.vdark_var.set(str(detector.get("dark_voltage_mv") or 0))
        self.band_min_var.set(str(detector.get("flat_band_min_nm", "")))
        self.band_max_var.set(str(detector.get("flat_band_max_nm", "")))
        if detector.get("responsivity_csv"):
            self._set_responsivity(os.path.join(self.base_dir, detector["responsivity_csv"]))
        self.detector_note.configure(text=detector.get("responsivity_note", ""))
        self._update_setup_note()

    def _update_setup_note(self):
        self.setup_note.configure(text=(
            f"No-film {self.v0_var.get()} mV · dark {self.vdark_var.get()} mV · "
            f"{self.band_min_var.get()}-{self.band_max_var.get()} nm · "
            + ("measurements imported from a file" if self._using_imported()
               else f"{len(self._measured_voltages())} built-in layer voltages")))

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
    def _imported_values(self):
        """Thicknesses that imported measurements sit at, so the graph can compare them."""
        if not self._using_imported():
            return []
        return [float(value) for value, _voltage in self.imported_measurements["points"]]

    def _sweep(self, peak):
        """Values plotted on the graph, always including any imported thicknesses."""
        imported = self._imported_values()
        if THICKNESS_MODES[self.thickness_mode_var.get()] == REFERENCE_MULTIPLIER:
            top = max([float(peak or 1.0)] + imported)
            values = [round(top * i / MIN_SWEEP_POINTS, 10) for i in range(1, MIN_SWEEP_POINTS + 1)]
        else:
            highest = max([float(peak or 1)] + imported)
            top = min(MAX_SWEEP_POINTS, max(MIN_SWEEP_POINTS, int(math.ceil(highest))))
            values = [float(i) for i in range(1, top + 1)]
        merged = sorted(set(values) | {round(value, 10) for value in imported})
        return merged[:MAX_SWEEP_POINTS * 2]

    def _entered_thickness(self):
        """The thickness the user typed, validated against the active thickness mode."""
        value = _required_float(self.value_var.get(), self.value_label.cget("text").split(". ", 1)[-1])
        if value <= 0:
            raise ValueError("Thickness must be greater than zero.")
        if THICKNESS_MODES[self.thickness_mode_var.get()] != REFERENCE_MULTIPLIER and not float(value).is_integer():
            raise ValueError("Layer counts must be whole numbers. For a fractional thickness switch the thickness "
                             "input to 'Reference-sample multiplier' under Advanced settings.")
        return value

    def build_config(self, sweep_peak=None):
        estimating = self.mode_var.get() == MODE_ESTIMATE
        detector = self.detectors[self._detector_key()]
        target_voltage = _required_float(self.value_var.get(), "Measured voltage") if estimating else None
        entered = None if estimating else self._entered_thickness()
        values = self._sweep(sweep_peak if sweep_peak is not None else entered)
        return BroadbandConfig(
            material=self.material_var.get(),
            values=values,
            no_film_voltage_mv=_required_float(self.v0_var.get(), "V0"),
            dark_voltage_mv=_optional_float(self.vdark_var.get(), "V_dark") or 0.0,
            thickness_mode=THICKNESS_MODES[self.thickness_mode_var.get()],
            reference_thickness_um=_optional_float(self.reference_thickness_var.get(), "x_ref"),
            layer_thickness_um=_optional_float(self.layer_thickness_var.get(), "Layer thickness"),
            layer_to_reference_ratio=_optional_float(self.ratio_var.get(), "Thickness ratio"),
            wavelength_min_nm=_optional_float(self.band_min_var.get(), "Window min"),
            wavelength_max_nm=_optional_float(self.band_max_var.get(), "Window max"),
            detector_name=detector.get("display_name", self._detector_key()),
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
            baseline_correction=BASELINE_MODES[self.baseline_var.get()],
            clamp_negative=self.clamp_var.get(),
            interface_correction=self.interface_var.get(),
            refractive_index=_optional_float(self.index_var.get(), "Refractive index"),
            measured_voltages_by_layer=(None if not self.compare_var.get() or self._using_imported()
                                        else self._measured_voltages()),
            measured_points=(self.imported_measurements["points"]
                             if self.compare_var.get() and self._using_imported() else None),
            measured_source=(self.imported_measurements["source"] if self._using_imported() else ""),
            inspect_value=entered if (entered is not None and entered in values) else None,
            target_voltage_mv=target_voltage,
        )

    def _estimated_peak(self, result):
        """Where the estimated thickness falls, so the graph sweep can cover it."""
        estimate = result.get("estimate")
        if estimate is None or estimate["status"] != "ok":
            return None
        if THICKNESS_MODES[self.thickness_mode_var.get()] == REFERENCE_MULTIPLIER:
            return estimate["thickness_scale"]
        return estimate["layer_equivalent"]

    def update_analysis(self):
        try:
            config = self.build_config()
            result = run_broadband_analysis(config, self.materials)
            peak = self._estimated_peak(result)
            # Estimate mode only learns the thickness from the first run; re-sweep around it.
            if peak and self._sweep(peak) != list(config.values):
                config = self.build_config(sweep_peak=peak)
                result = run_broadband_analysis(config, self.materials)
            self.result = result
        except Exception as exc:
            messagebox.showerror("Analysis error", str(exc), parent=self)
            return
        self._update_setup_note()
        self.baseline_note.configure(text=self.result["corrections"]["baseline_reason"])
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

    def _estimate_x(self):
        """The estimated thickness expressed on the current graph x-axis."""
        estimate = self.result.get("estimate")
        if estimate is None:
            return None
        mode = self.result["thickness"]["mode"]
        if mode == LAYERS_PHYSICAL:
            return estimate["thickness_um"]
        if mode == LAYERS_RATIO:
            return estimate["layer_equivalent"]
        return estimate["thickness_scale"]

    # ------------------------------------------------------------- reporting
    def _detector_label(self):
        return self.detectors[self._detector_key()].get("display_name", self._detector_key())

    def _thickness_words(self, row):
        if row["thickness_um"] is not None:
            return f"{row['layer_count']} layer(s) = {row['thickness_um']:g} µm"
        if row["layer_count"] is not None:
            return f"{row['layer_count']} layer(s)"
        return f"{row['thickness_scale']:g}× the reference sample"

    def _prediction_text(self):
        """Big answer plus one line of context, for predict mode."""
        result = self.result
        index = result["inspected_index"]
        if index is None:
            return "", ""
        row, corrected = result["rows"][index], result["corrections"]["any_enabled"]
        voltage = row["predicted_voltage_mv"] if corrected else row["raw_predicted_voltage_mv"]
        parts = [f"{self._thickness_words(row)} of {result['material']['name']}", self._detector_label()]
        if row["measured_voltage_mv"] is not None:
            direction = "above" if row["error_mv"] > 0 else "below"
            parts.append(f"measured {row['measured_voltage_mv']:.1f} mV — the model sits "
                         f"{abs(row['percent_error']):.0f}% {direction} it")
        if corrected:
            parts.append(f"raw model {row['raw_predicted_voltage_mv']:.1f} mV")
        return f"{voltage:.1f} mV predicted", "   ·   ".join(parts)

    def _estimate_text(self):
        """Big answer plus one line of context, for estimate mode."""
        estimate = self.result["estimate"]
        if estimate["thickness_um"] is not None:
            headline = f"{estimate['thickness_um']:.1f} µm"
        elif estimate["layer_equivalent"] is not None:
            headline = f"{estimate['layer_equivalent']:.2f} layers"
        else:
            headline = f"{estimate['thickness_scale']:.3f}× the reference sample"
        parts = [f"from {estimate['target_voltage_mv']:.1f} mV measured", self.result["material"]["name"],
                 self._detector_label()]
        solved = estimate["status"] == "ok"
        if solved and estimate["thickness_um"] is None and estimate["layer_equivalent"] is not None:
            parts.append(f"x/x_ref = {estimate['thickness_scale']:.3f}")
        if solved and estimate["partial_coverage"]:
            low, high = estimate["thickness_scale_bounds"]
            parts.append(f"x/x_ref could be {low:.2f}–{'unbounded' if high is None else f'{high:.2f}'} "
                         f"given the spectral coverage")
        if estimate["status"] == "exceeds_max_scale":
            headline = "No thickness fits"
            parts.insert(0, "this voltage is darker than the model can ever reach, so the film is thicker than "
                            "this measurement can resolve (a LOWER limit only)")
        elif estimate["status"] == "at_or_above_no_film":
            headline = "0 µm (upper limit)"
            parts.insert(0, "this voltage is at or above the no-film prediction, so the model needs no film at all")
        return headline, "   ·   ".join(parts)

    def _status_line(self):
        """One short line; the full wording lives under 'Numbers'."""
        result = self.result
        corrections = result["corrections"]
        problems = [w for w in result["warnings"] if not w.startswith("BASELINE CORRECTED AUTOMATICALLY")]
        if problems:
            first = problems[0].split(". ")[0].strip()
            if len(first) > 150:
                first = first[:147] + "…"
            extra = f"   (+{len(problems) - 1} more)" if len(problems) > 1 else ""
            return (f"⚠  {first}{extra}   ·   open 'Numbers' for the full wording",
                    "#7f1d1d" if result["coverage"]["partial"] else "#78350f")
        if corrections["baseline_automatic"]:
            return (f"Baseline corrected automatically: this spectrum's transparent region sat "
                    f"{corrections['baseline_offset']:+.4f} absorbance below zero, so that offset was subtracted "
                    f"before the model ran.", "#1e3a8a")
        return (f"Coverage {result['coverage']['fraction']:.0%} of the detector band   ·   no warnings",
                "#14532d")

    def _forward_model_lines(self, scale=None, t_eff=None, voltage=None):
        """The forward model in words, with this run's numbers when there are any."""
        result = self.result
        corrections = result["corrections"]
        v0, vd = result["voltages"]["no_film_voltage_mv"], result["voltages"]["dark_voltage_mv"]
        low, high = result["coverage"]["weight_range_nm"]
        ratio = "x / x_ref" if scale is None else f"{scale:g}"
        steps = []
        if corrections["baseline_enabled"]:
            why = ("it sat below zero, which would make the film transmit more than 100 %"
                   if corrections["baseline_automatic"] else "the correction is set to always subtract it")
            steps.append([f"{result['material']['name']}'s flat baseline of {corrections['baseline_offset']:+.4f} "
                          f"absorbance is subtracted first, because {why}.",
                          "   The absorption peaks themselves are untouched."])
        steps.append([f"The measured absorbance A(λ) is scaled to the thickness:  A(λ) × {ratio}",
                      "   (x / x_ref = this film divided by the film the spectrum was measured on)."])
        steps.append([f"Each wavelength becomes a transmission:  T(λ) = 10^(−A(λ) × {ratio})."])
        weighted = [f"T(λ) is averaged over {low:.0f}–{high:.0f} nm, weighted by source × detector "
                    f"response" + (f":  T_eff = {t_eff:.4f}." if t_eff is not None else "  →  T_eff.")]
        if corrections["interface_enabled"]:
            weighted.append(f"   Fresnel reflection at each layer's two surfaces is included: "
                            f"×{corrections['interface_per_layer']:.4f} per layer.")
        steps.append(weighted)
        last = "Voltage = V_dark + (V0 − V_dark) × T_eff"
        if voltage is not None:
            last += f" = {vd:g} + ({v0:g} − {vd:g}) × {t_eff:.4f} = {voltage:.1f} mV."
        else:
            last += f",  with V0 = {v0:g} mV and V_dark = {vd:g} mV."
        steps.append([last])
        lines = ["The forward model — nothing in it is fitted to the measured voltages:"]
        for number, step in enumerate(steps, start=1):
            lines.append(f"{number}. {step[0]}")
            lines.extend(step[1:])
        return lines

    def _explanation_text(self):
        """Plain-language walk-through of the active calculation."""
        result = self.result
        v0, vd = result["voltages"]["no_film_voltage_mv"], result["voltages"]["dark_voltage_mv"]
        estimate = result.get("estimate")
        if estimate is None:
            index = result["inspected_index"]
            if index is None:
                return ""
            row = result["rows"][index]
            lines = self._forward_model_lines(row["thickness_scale"], row["predicted_transmission"],
                                              row["predicted_voltage_mv"])
            if row["measured_voltage_mv"] is not None:
                lines.append(f"{sum(1 for x in lines if x[:1].isdigit()) + 1}. Measured "
                             f"{row['measured_voltage_mv']:.1f} mV, so the model is off by "
                             f"{row['error_mv']:+.1f} mV ({row['percent_error']:+.1f}%).")
            return "\n".join(lines)

        solved = estimate["status"] == "ok"
        lines = ["Estimating a thickness = the same forward model, run backwards.",
                 f"1. The measurement becomes a transmission:  T_eff = (V \u2212 V_dark) / (V0 \u2212 V_dark)"
                 f" = ({estimate['target_voltage_mv']:.1f} \u2212 {vd:g}) / ({v0:g} \u2212 {vd:g}) = "
                 f"{estimate['target_transmission']:.4f}",
                 "2. Trial thicknesses are run through the forward model below."]
        if solved:
            lines += [f"3. The thickness is bracketed and halved {estimate['iterations']} times until the predicted "
                      f"voltage matches:  x / x_ref = {estimate['thickness_scale']:.4f}",
                      f"4. Check: that thickness predicts {estimate['back_predicted_voltage_mv']:.3f} mV, "
                      f"{estimate['residual_mv']:+.1e} mV away from the measurement.", ""]
            lines += self._forward_model_lines(estimate["thickness_scale"], estimate["effective_transmission"])
        else:
            # No thickness reproduces the measurement, so there are no numbers to quote.
            lines += ["3. No thickness reproduces that voltage.",
                      f"   {estimate['status_note']}", ""]
            lines += self._forward_model_lines()
        return "\n".join(lines)

    def _update_status(self):
        headline, subline = self._estimate_text() if self.result.get("estimate") else self._prediction_text()
        self.headline_label.configure(text=headline)
        self.subline_label.configure(text=subline)
        text, color = self._status_line()
        self.status_label.configure(text=text, fg_color=color)
        self.explain_label.configure(text=self._explanation_text())

    def _update_results(self):
        result = self.result
        fmt = lambda v, spec: "—" if v is None else format(v, spec)
        lines = [f"WARNING: {warning}" for warning in result["warnings"]]
        if lines:
            lines.append("")
        lines += [f"Material:  {result['material']['name']} — {result['material']['note']}",
                  f"Thickness: {result['thickness']['description']}",
                 f"Detector:  {result['weighting']['detector_description']}",
                 f"Source:    {result['weighting']['source_description']}",
                 f"V0 = {result['voltages']['no_film_voltage_mv']:g} mV, V_dark = {result['voltages']['dark_voltage_mv']:g} mV;  "
                 f"Vpred = V_dark + (V0 - V_dark)·T_eff"]
        estimate = result.get("estimate")
        if estimate:
            lines += ["",
                      f"Thickness estimate: measured {estimate['target_voltage_mv']:g} mV -> T_eff = "
                      f"{estimate['target_transmission']:.4f} -> x/x_ref = {estimate['thickness_scale']:.4f} "
                      f"({fmt(estimate['layer_equivalent'], '.3f')} layers, {fmt(estimate['thickness_um'], '.2f')} µm)",
                      f"  bisection of the forward model: {estimate['iterations']} iterations, "
                      f"converged={estimate['converged']}, back-predicted "
                      f"{estimate['back_predicted_voltage_mv']:.4f} mV (residual {estimate['residual_mv']:+.2e} mV)"]
        lines += ["",
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
    def _uncorrected_label(self):
        """Name the reference curve after whatever the corrections actually changed."""
        corrections = self.result["corrections"]
        if corrections["baseline_enabled"]:
            return "Before baseline correction"
        return "Predicted (raw model)"

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

    def _mark_estimate(self, show_voltage=True):
        """Cross-hairs at the measured voltage and the thickness solved from it."""
        estimate = self.result.get("estimate")
        x = self._estimate_x()
        if estimate is None or x is None or estimate["status"] != "ok":
            return
        if show_voltage:
            self.ax.axhline(estimate["target_voltage_mv"], color="#fb7185", linestyle="--", linewidth=1.2,
                            label=f"Measured {estimate['target_voltage_mv']:g} mV")
        self.ax.axvline(x, color="#fbbf24", linestyle="--", linewidth=1.4, label=f"Estimated thickness ({x:.3g})")

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
        elif graph == GRAPH_VOLTAGE:
            self._style_axis(graph, xlabel, "Detector voltage (mV)")
            self.ax.plot(xs, [r["raw_predicted_voltage_mv"] for r in rows], marker="o", color="#34d399",
                         label=self._uncorrected_label())
            if corrected:
                self.ax.plot(xs, [r["predicted_voltage_mv"] for r in rows], marker="s", color="#38bdf8",
                             label="Predicted (with corrections)")
            if result["coverage"]["partial"]:
                self.ax.fill_between(xs, [r["predicted_voltage_bounds_mv"][0] for r in rows],
                                     [r["predicted_voltage_bounds_mv"][1] for r in rows], color="#f87171", alpha=0.15,
                                     label="Bounds (uncovered weight T=0…1)")
            if compared:
                self.ax.scatter([x for x, _ in compared], [r["measured_voltage_mv"] for _, r in compared],
                                color="#fb7185", s=70, zorder=3, label="Measured")
            self._mark_estimate()
            self._legend()
        elif graph == GRAPH_PARITY:
            self._style_axis(graph, "Measured voltage (mV)", "Predicted voltage (mV)")
            if compared:
                measured = [r["measured_voltage_mv"] for _, r in compared]
                keys = [("raw_predicted_voltage_mv", self._uncorrected_label(), "#34d399")]
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
            self.ax.plot(xs, [r["raw_effective_transmission"] for r in rows], marker="o", color="#34d399",
                         label=self._uncorrected_label())
            if corrected:
                self.ax.plot(xs, [r["predicted_transmission"] for r in rows], marker="s", color="#38bdf8", label="Corrected")
            if compared:
                self.ax.scatter([x for x, _ in compared], [r["measured_transmission"] for _, r in compared],
                                color="#fb7185", s=60, zorder=3, label="Measured")
            self._mark_estimate(show_voltage=False)
            self._legend()
        else:
            self._style_axis(graph, xlabel, "ln(1/T_eff)")
            self.ax.plot(xs, [np.nan if r["predicted_attenuation_natural"] is None else r["predicted_attenuation_natural"]
                              for r in rows], marker="o", color="#a78bfa", label="Predicted (corrections as enabled)")
            if compared:
                self.ax.scatter([x for x, _ in compared],
                                [np.nan if r["measured_attenuation_natural"] is None else r["measured_attenuation_natural"]
                                 for _, r in compared], color="#fb7185", s=60, zorder=3, label="Measured")
            self._mark_estimate(show_voltage=False)
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
        default_name = f"broadband_{self.material_var.get()}_{self._detector_key()}.xlsx"
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
