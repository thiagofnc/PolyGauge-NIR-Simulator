"""CustomTkinter window for broadband detector analysis and export."""

from __future__ import annotations

import json
import math
import os
from tkinter import filedialog, messagebox

import customtkinter as ctk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from BroadbandAnalysis import (analyze_layers, estimate_baseline_offset,
                               fresnel_layer_transmission)
from BroadbandExcel import export_broadband_workbook
from MeasuredData import discover_measured_samples, load_log_spectrum, wavenumber_to_nm


GRAPH_ABSORBANCE = "Absorbance Spectrum (raw vs corrected)"
GRAPH_TRANSMISSION = "Transmission Spectrum per Layer Count"
GRAPH_WEIGHTING = "Detector Weighting W(λ)"
GRAPH_CONTRIBUTION = "Wavelengths Contributing to Detector"
GRAPH_VOLTAGE = "Predicted Voltage vs Layer Count"
GRAPH_COMPARE = "Measured vs Predicted Voltage"
GRAPH_PARITY = "Measured vs Predicted (Parity Plot)"
GRAPH_ERROR = "Prediction Error vs Layer Count"
GRAPH_EFFECTIVE_T = "Effective Transmission vs Layer Count"
GRAPH_ATTENUATION = "Effective Attenuation vs Layer Count"
GRAPH_TYPES = [
    GRAPH_ABSORBANCE, GRAPH_TRANSMISSION, GRAPH_WEIGHTING, GRAPH_CONTRIBUTION,
    GRAPH_VOLTAGE, GRAPH_COMPARE, GRAPH_PARITY, GRAPH_ERROR, GRAPH_EFFECTIVE_T,
    GRAPH_ATTENUATION,
]

BASELINE_NONE = "None (use spectrum as recorded)"
BASELINE_PERCENTILE = "Subtract in-band baseline (5th percentile)"

# Typical mid-IR refractive indices used for per-layer Fresnel reflection loss.
DEFAULT_REFRACTIVE_INDEX = {"Nylon": 1.53, "EVOH": 1.52, "PE": 1.51}
LAYER_COLORS = ["#38bdf8", "#34d399", "#fbbf24", "#fb7185", "#a78bfa", "#f472b6", "#2dd4bf", "#fb923c"]


def _sorted_spectrum(wavelength_nm, values):
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(wavelength_nm) & np.isfinite(values) & (wavelength_nm > 0)
    order = np.argsort(wavelength_nm[valid])
    return wavelength_nm[valid][order], values[valid][order]


def load_analysis_materials(base_dir):
    """Load named spectra plus every absorbance log discovered by the app."""
    materials = {}
    known_logs = {
        "Nylon": os.path.join(base_dir, "logs_full_range", "nylon_reynolds_abs.txt"),
        "PE": os.path.join(base_dir, "logs_full_range", "PE_reynolds_abs.txt"),
    }
    for name, path in known_logs.items():
        if os.path.exists(path):
            wavenumber, values, header = load_log_spectrum(path)
            wavelength, values = _sorted_spectrum(wavenumber_to_nm(wavenumber), values)
            materials[name] = {
                "wavelength_nm": wavelength, "absorbance": values, "mode": "base10",
                "source": os.path.relpath(path, base_dir),
                "note": f"Measured FTIR absorbance ({header.get('YUNITS', 'Abs')}); treated as one reference layer.",
            }

    evoh_path = os.path.join(base_dir, "EVOHTransmissionPercentvsWavelength_in_nm_2000nm_to_5000nm.csv")
    if os.path.exists(evoh_path):
        frame = pd.read_csv(evoh_path, skipinitialspace=True)
        wavelength, transmission = _sorted_spectrum(frame.iloc[:, 0], frame.iloc[:, 1])
        with np.errstate(divide="ignore", invalid="ignore"):
            absorbance = -np.log10(np.clip(transmission, np.finfo(float).tiny, None))
        materials["EVOH"] = {
            "wavelength_nm": wavelength, "absorbance": absorbance, "mode": "base10",
            "source": os.path.basename(evoh_path),
            "note": "Digitized EVOH transmission converted with A=-log10(T); treated as one reference layer.",
        }

    samples, _headers = discover_measured_samples(
        os.path.join(base_dir, "logs_full_range"), os.path.join(base_dir, "sample_references.xlsx")
    )
    for sample in samples:
        path = sample["paths"].get("absorbance")
        if not path:
            continue
        display = sample["label"]
        if any(os.path.normcase(path) == os.path.normcase(item["source"] if os.path.isabs(item["source"]) else os.path.join(base_dir, item["source"]))
               for item in materials.values()):
            continue
        wavenumber, values, header = load_log_spectrum(path)
        wavelength, values = _sorted_spectrum(wavenumber_to_nm(wavenumber), values)
        materials[f"Measured: {display}"] = {
            "wavelength_nm": wavelength, "absorbance": values, "mode": "base10",
            "source": os.path.relpath(path, base_dir),
            "note": f"Loaded measured absorbance ({header.get('YUNITS', 'Abs')}); reference thickness must be supplied by the user.",
        }
    return materials


def _blackbody(wavelength_nm, temperature_k):
    h, c, k = 6.62607015e-34, 299792458.0, 1.380649e-23
    wl_m = np.asarray(wavelength_nm, dtype=float) * 1e-9
    exponent = np.clip(h * c / (wl_m * k * temperature_k), 0, 700)
    intensity = (2 * h * c ** 2 / wl_m ** 5) / np.expm1(exponent)
    maximum = np.nanmax(intensity)
    return intensity / maximum if maximum > 0 else np.zeros_like(wl_m)


class BroadbandAnalysisWindow(ctk.CTkToplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        self.title("Broadband Spectral Detector Analysis")
        self.geometry("1500x900")
        self.minsize(1120, 700)
        self.materials = load_analysis_materials(self.base_dir)
        with open(os.path.join(self.base_dir, "experimental_detector_data.json"), encoding="utf-8") as handle:
            self.detectors = json.load(handle)
        if not self.materials:
            messagebox.showerror("No spectra", "No absorbance spectra could be loaded.", parent=self)
            self.destroy()
            return
        self.analysis = None
        self._build_ui()
        self._material_changed()
        self._detector_changed()
        self.update_analysis()

    def _build_ui(self):
        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
        controls = ctk.CTkScrollableFrame(self, width=340, label_text="Analysis Settings")
        controls.grid(row=0, column=0, sticky="nsew", padx=(10, 5), pady=10)
        view = ctk.CTkFrame(self)
        view.grid(row=0, column=1, sticky="nsew", padx=(5, 10), pady=10)

        def heading(text):
            ctk.CTkLabel(controls, text=text, anchor="w", font=("Arial", 15, "bold"),
                         text_color="#7dd3fc").pack(fill="x", padx=8, pady=(16, 0))

        def option(label, variable, values, command=None):
            ctk.CTkLabel(controls, text=label, anchor="w", font=("Arial", 13, "bold")).pack(fill="x", padx=8, pady=(10, 2))
            widget = ctk.CTkOptionMenu(controls, variable=variable, values=values, command=command)
            widget.pack(fill="x", padx=8)
            return widget

        def entry(label, variable):
            ctk.CTkLabel(controls, text=label, anchor="w").pack(fill="x", padx=8, pady=(9, 1))
            ctk.CTkEntry(controls, textvariable=variable).pack(fill="x", padx=8)

        heading("1. Material")
        self.material_var = ctk.StringVar(value=next(iter(self.materials)))
        option("Material spectrum", self.material_var, list(self.materials), self._material_changed)

        heading("2. Detector & Layers")
        self.detector_var = ctk.StringVar(value=next(iter(self.detectors)))
        option("Detector", self.detector_var, list(self.detectors), lambda _value: self._detector_changed())
        self.scale_mode_var = ctk.StringVar(value="Layer count")
        option("Simulation input", self.scale_mode_var, ["Layer count", "Thickness multiplier"])
        range_frame = ctk.CTkFrame(controls, fg_color="transparent")
        range_frame.pack(fill="x", padx=8, pady=(10, 0))
        self.start_var, self.end_var, self.step_var = (ctk.StringVar(value=value) for value in ("1", "4", "1"))
        for label, variable in [("Start", self.start_var), ("End", self.end_var), ("Step", self.step_var)]:
            part = ctk.CTkFrame(range_frame, fg_color="transparent")
            part.pack(side="left", expand=True, fill="x", padx=2)
            ctk.CTkLabel(part, text=label).pack()
            ctk.CTkEntry(part, textvariable=variable, width=75).pack(fill="x")
        self.selected_scale_var = ctk.StringVar(value="3")
        entry("Layers to inspect (headline + spectral graphs)", self.selected_scale_var)
        self.no_film_var = ctk.StringVar(value="498")
        self.min_wl_var = ctk.StringVar(value="2000")
        self.max_wl_var = ctk.StringVar(value="12000")
        entry("No-film voltage V0 (mV)", self.no_film_var)
        entry("Detector minimum wavelength (nm)", self.min_wl_var)
        entry("Detector maximum wavelength (nm)", self.max_wl_var)
        self.weighting_var = ctk.StringVar(value="Flat wavelength range")
        option("Detector weighting", self.weighting_var,
               ["Flat wavelength range", "Detector responsivity", "Source × detector responsivity"])
        self.temp_var = ctk.StringVar(value="3000")
        entry("Source temperature (K, source weighting only)", self.temp_var)

        heading("3. Physics Model")
        self.mode_var = ctk.StringVar(value="Base-10 absorbance")
        option("Absorbance convention", self.mode_var,
               ["Base-10 absorbance", "Natural-log / Napierian attenuation"])
        self.baseline_var = ctk.StringVar(value=BASELINE_PERCENTILE)
        option("Baseline correction", self.baseline_var, [BASELINE_PERCENTILE, BASELINE_NONE])
        self.clamp_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(controls, text="Clamp negative absorbance to zero", variable=self.clamp_var).pack(anchor="w", padx=8, pady=(12, 4))
        self.fresnel_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(controls, text="Include Fresnel reflection loss per layer", variable=self.fresnel_var).pack(anchor="w", padx=8, pady=4)
        self.refractive_index_var = ctk.StringVar(value="1.5")
        entry("Film refractive index n", self.refractive_index_var)
        self.reference_thickness_var = ctk.StringVar(value="1 reference layer")
        entry("Reference thickness / description", self.reference_thickness_var)
        self.compare_var = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(controls, text="Compare with measured detector voltages", variable=self.compare_var).pack(anchor="w", padx=8, pady=(10, 4))

        self.warning_label = ctk.CTkLabel(controls, text="", text_color="#ffcc66", wraplength=310, justify="left")
        self.warning_label.pack(fill="x", padx=8, pady=6)
        self.note_label = ctk.CTkLabel(controls, text="", text_color="#a9c8e8", wraplength=310, justify="left")
        self.note_label.pack(fill="x", padx=8, pady=4)
        ctk.CTkButton(controls, text="UPDATE ANALYSIS", command=self.update_analysis,
                      fg_color="#00aa00", hover_color="#008800", font=("Arial", 14, "bold")).pack(fill="x", padx=8, pady=(10, 5))
        ctk.CTkButton(controls, text="EXPORT EXCEL", command=self.export_excel,
                      fg_color="#1f6f43", hover_color="#278b55", font=("Arial", 14, "bold")).pack(fill="x", padx=8, pady=5)

        view.grid_rowconfigure(2, weight=1)
        view.grid_columnconfigure(0, weight=1)
        self.headline_label = ctk.CTkLabel(view, text="", font=("Arial", 20, "bold"), anchor="w", justify="left")
        self.headline_label.grid(row=0, column=0, sticky="ew", padx=12, pady=(10, 0))
        graph_bar = ctk.CTkFrame(view, fg_color="transparent")
        graph_bar.grid(row=1, column=0, sticky="ew", padx=8, pady=(6, 0))
        ctk.CTkLabel(graph_bar, text="Graph:", font=("Arial", 13, "bold")).pack(side="left", padx=(4, 6))
        self.graph_var = ctk.StringVar(value=GRAPH_COMPARE)
        ctk.CTkOptionMenu(graph_bar, variable=self.graph_var, values=GRAPH_TYPES, width=340,
                          command=lambda _value: self.draw_graph()).pack(side="left")
        self.fig, self.ax = plt.subplots(figsize=(10, 6), facecolor="#2b2b2b")
        self.canvas = FigureCanvasTkAgg(self.fig, master=view)
        self.canvas.get_tk_widget().grid(row=2, column=0, sticky="nsew", padx=8, pady=8)
        self.results_text = ctk.CTkTextbox(view, height=160, font=("Consolas", 12))
        self.results_text.grid(row=3, column=0, sticky="ew", padx=8, pady=(0, 8))

    def _material_changed(self, _value=None):
        name = self.material_var.get()
        material = self.materials[name]
        self.mode_var.set("Base-10 absorbance" if material["mode"] == "base10" else "Natural-log / Napierian attenuation")
        self.refractive_index_var.set(str(DEFAULT_REFRACTIVE_INDEX.get(name, 1.5)))

    def _detector_changed(self):
        detector = self.detectors[self.detector_var.get()]
        self.no_film_var.set(str(detector["no_film_voltage_mv"]))
        self.min_wl_var.set(str(detector["wavelength_min_nm"]))
        self.max_wl_var.set(str(detector["wavelength_max_nm"]))

    def _scales(self):
        start, end, step = float(self.start_var.get()), float(self.end_var.get()), float(self.step_var.get())
        if step <= 0 or end < start:
            raise ValueError("Simulation range requires end ≥ start and step > 0.")
        count = int(math.floor((end - start) / step + 1e-10)) + 1
        if count > 1000:
            raise ValueError("Simulation range is limited to 1000 points.")
        return [start + i * step for i in range(count)]

    def _weight_curves(self, wavelength, lower, upper):
        mode = self.weighting_var.get()
        response = None
        source = None
        if mode != "Flat wavelength range":
            response = (np.array([lower, upper]), np.ones(2))
        if mode == "Source × detector responsivity":
            source = (wavelength, _blackbody(wavelength, float(self.temp_var.get())))
        return source, response

    def update_analysis(self):
        try:
            material_name = self.material_var.get()
            material = self.materials[material_name]
            detector = self.detectors[self.detector_var.get()]
            lower, upper = float(self.min_wl_var.get()), float(self.max_wl_var.get())
            wavelength, raw_absorbance = material["wavelength_nm"], material["absorbance"]
            baseline_offset = 0.0
            if self.baseline_var.get() == BASELINE_PERCENTILE:
                baseline_offset = estimate_baseline_offset(wavelength, raw_absorbance, lower, upper)
            interface = (fresnel_layer_transmission(float(self.refractive_index_var.get()))
                         if self.fresnel_var.get() else 1.0)
            source, response = self._weight_curves(wavelength, lower, upper)
            measured = {}
            if self.compare_var.get():
                raw = detector.get("materials", {}).get(material_name, {})
                measured = {float(key): float(value) for key, value in raw.items()}
            absorbance_mode = "base10" if self.mode_var.get().startswith("Base-10") else "natural"
            self.analysis = analyze_layers(
                wavelength, raw_absorbance - baseline_offset, float(self.no_film_var.get()),
                self._scales(), measured, float(self.selected_scale_var.get()), source=source,
                responsivity=response, wavelength_min=lower, wavelength_max=upper,
                absorbance_mode=absorbance_mode, clamp_negative=self.clamp_var.get(),
                layer_interface_transmission=interface,
            )
            self.analysis["baseline_offset"] = baseline_offset
            self.analysis["layer_interface_transmission"] = interface
            spectral = self.analysis["spectral"]
            spectral["raw_absorbance"] = np.interp(spectral["wavelength_nm"], wavelength, raw_absorbance)
            actual_min, actual_max = spectral["integration_range_nm"]
            warnings = []
            if np.any(spectral["absorbance"] < 0) and not self.clamp_var.get():
                warnings.append("Negative absorbance is present in band (baseline/noise); predicted T may exceed 1.")
            if self.weighting_var.get() == "Flat wavelength range":
                warnings.append("Flat in-band weighting is an approximation.")
            else:
                warnings.append("Detector response is a rectangular approximation; no measured response curve is bundled.")
            if actual_min > lower or actual_max < upper:
                warnings.append(f"Spectrum only covers {actual_min:.0f}–{actual_max:.0f} nm of the detector band; "
                                "absorption outside that range is not counted.")
            self.warning_label.configure(text="\n".join(warnings))
            self.note_label.configure(text=material["note"] + "\n" + detector["responsivity_note"])
            self._update_results()
            self.draw_graph()
        except Exception as exc:
            messagebox.showerror("Analysis error", str(exc), parent=self)

    def _selected_row(self):
        selected = self.analysis["selected_scale"]
        return next(row for row in self.analysis["rows"] if row["thickness_scale"] == selected)

    def _scale_label(self):
        return "Layer count" if self.scale_mode_var.get() == "Layer count" else "Thickness multiplier"

    def _update_results(self):
        row = self._selected_row()
        if self.scale_mode_var.get() == "Layer count":
            amount = f"{row['thickness_scale']:g} layer(s)"
        else:
            amount = f"×{row['thickness_scale']:g} thickness"
        headline = f"{self.material_var.get()}, {amount}:   Predicted {row['predicted_voltage_mv']:.1f} mV"
        if row["measured_voltage_mv"] is not None:
            percent = "—" if row["percent_error"] is None else f"{row['percent_error']:+.1f}%"
            headline += (f"   |   Measured {row['measured_voltage_mv']:.1f} mV"
                         f"   |   Error {row['error_mv']:+.1f} mV / {percent}")
        self.headline_label.configure(text=headline)

        self.results_text.configure(state="normal")
        self.results_text.delete("1.0", "end")
        metrics = self.analysis["metrics"]
        lines = ["Layer   Measured mV   Predicted mV   Spectral T   Interface T   Total T   Error mV   Error %"]
        for row in self.analysis["rows"]:
            measured = "—" if row["measured_voltage_mv"] is None else f"{row['measured_voltage_mv']:.1f}"
            error = "—" if row["error_mv"] is None else f"{row['error_mv']:+.1f}"
            percent = "—" if row["percent_error"] is None else f"{row['percent_error']:+.1f}"
            lines.append(f"{row['layer_count']!s:<7} {measured:>11}   {row['predicted_voltage_mv']:>12.1f}   "
                         f"{row['spectral_transmission']:>10.4f}   {row['interface_transmission']:>11.4f}   "
                         f"{row['predicted_transmission']:>7.4f}   {error:>8}   {percent:>7}")
        if metrics["mae_mv"] is not None:
            r2 = "—" if metrics["r_squared"] is None else f"{metrics['r_squared']:.4f}"
            mape = "—" if metrics["mape_percent"] is None else f"{metrics['mape_percent']:.1f}%"
            lines.append(f"\nModel agreement: MAE={metrics['mae_mv']:.2f} mV   RMSE={metrics['rmse_mv']:.2f} mV   "
                         f"MAPE={mape}   R²={r2}")
        lines.append(f"Baseline offset removed: {self.analysis['baseline_offset']:+.4f} A   "
                     f"Interface transmission per layer: {self.analysis['layer_interface_transmission']:.4f}")
        self.results_text.insert("1.0", "\n".join(lines))
        self.results_text.configure(state="disabled")

    def _style_axis(self, title, xlabel, ylabel):
        self.ax.clear()
        self.ax.set_facecolor("#2b2b2b")
        self.ax.tick_params(colors="white")
        for spine in self.ax.spines.values():
            spine.set_color("gray")
        self.ax.set_title(title, color="white")
        self.ax.set_xlabel(xlabel, color="white")
        self.ax.set_ylabel(ylabel, color="white")
        self.ax.grid(True, alpha=0.2)

    def _legend(self, title=None):
        self.ax.legend(title=title, facecolor="#333333", edgecolor="gray", labelcolor="white")

    def draw_graph(self):
        if self.analysis is None:
            return
        graph = self.graph_var.get()
        spectral = self.analysis["spectral"]
        rows = self.analysis["rows"]
        scales = [r["thickness_scale"] for r in rows]
        compared = [r for r in rows if r["measured_voltage_mv"] is not None]
        compared_scales = [r["thickness_scale"] for r in compared]
        x = spectral["wavelength_nm"]
        scale_label = self._scale_label()
        selected = self.analysis["selected_scale"]

        if graph == GRAPH_ABSORBANCE:
            self._style_axis(f"{self.material_var.get()} absorbance (one reference layer)", "Wavelength (nm)", "Absorbance")
            if self.analysis["baseline_offset"]:
                self.ax.plot(x, spectral["raw_absorbance"], color="#94a3b8", linewidth=1, label="Raw spectrum")
            self.ax.plot(x, spectral["absorbance"], color="#38bdf8", label="Used by model")
            self.ax.axhline(0, color="gray", linewidth=0.8)
            self._legend()
        elif graph == GRAPH_TRANSMISSION:
            self._style_axis("Transmission T(λ) through the film stack", "Wavelength (nm)", "Transmission")
            spectra = self.analysis["spectra_by_scale"]
            for i, (scale, result) in enumerate(spectra.items()):
                if len(spectra) > len(LAYER_COLORS) and scale != selected:
                    continue
                self.ax.plot(x, result["transmission"] * result["interface_transmission"],
                             color=LAYER_COLORS[i % len(LAYER_COLORS)], label=f"{scale:g}",
                             linewidth=2.2 if scale == selected else 1.2)
            self._legend(scale_label)
        elif graph == GRAPH_WEIGHTING:
            self._style_axis("Detector weighting", "Wavelength (nm)", "Normalized spectral value")
            if spectral["source"] is not None:
                self.ax.plot(x, spectral["source"], label="Source S(λ)", color="#ffcc00")
            if spectral["responsivity"] is not None:
                self.ax.plot(x, spectral["responsivity"], label="Detector R(λ)", color="#38bdf8")
            peak = np.max(np.abs(spectral["weight"]))
            normalized = spectral["weight"] / peak if peak else spectral["weight"]
            self.ax.plot(x, normalized, label="Normalized W(λ)", color="#f472b6", linewidth=2)
            self._legend()
        elif graph == GRAPH_CONTRIBUTION:
            peak = np.max(np.abs(spectral["weight"])) or 1.0
            self._style_axis(f"What the detector sees at {selected:g} ({scale_label.lower()})",
                             "Wavelength (nm)", "Normalized contribution")
            self.ax.fill_between(x, spectral["weight"] / peak, color="#38bdf8", alpha=0.2, label="No film: W(λ)")
            self.ax.fill_between(x, spectral["weighted_transmitted"] * spectral["interface_transmission"] / peak,
                                 color="#fb7185", alpha=0.6, label="With film: W(λ)·T(λ)")
            self.ax.text(0.01, 0.03, f"Surviving signal: {spectral['spectral_transmission']:.1%} spectral × "
                         f"{spectral['interface_transmission']:.1%} interface = {spectral['effective_transmission']:.1%}",
                         transform=self.ax.transAxes, color="white")
            self._legend()
        elif graph == GRAPH_VOLTAGE:
            self._style_axis(graph, scale_label, "Detector voltage (mV)")
            self.ax.plot(scales, [r["predicted_voltage_mv"] for r in rows], marker="o", color="#34d399")
        elif graph == GRAPH_COMPARE:
            self._style_axis(graph, scale_label, "Detector voltage (mV)")
            self.ax.plot(scales, [r["predicted_voltage_mv"] for r in rows], marker="o",
                         label="Predicted (spectral model)", color="#34d399")
            if compared:
                self.ax.scatter(compared_scales, [r["measured_voltage_mv"] for r in compared],
                                label="Measured", color="#fb7185", s=60, zorder=3)
            self._legend()
        elif graph == GRAPH_PARITY:
            self._style_axis(graph, "Measured voltage (mV)", "Predicted voltage (mV)")
            if compared:
                measured = [r["measured_voltage_mv"] for r in compared]
                predicted = [r["predicted_voltage_mv"] for r in compared]
                low, high = min(measured + predicted) * 0.95, max(measured + predicted) * 1.05
                self.ax.plot([low, high], [low, high], color="gray", linestyle="--", label="Perfect agreement")
                self.ax.scatter(measured, predicted, color="#fbbf24", s=60, zorder=3, label=scale_label)
                for r in compared:
                    self.ax.annotate(f"{r['thickness_scale']:g}", (r["measured_voltage_mv"], r["predicted_voltage_mv"]),
                                     textcoords="offset points", xytext=(6, 4), color="white")
                self._legend()
            else:
                self.ax.text(0.5, 0.5, "No measured data for this material/detector", color="white",
                             ha="center", transform=self.ax.transAxes)
        elif graph == GRAPH_ERROR:
            self._style_axis(graph, scale_label, "Prediction error (%)  [+ = over-predicts]")
            errors = [r["percent_error"] for r in compared if r["percent_error"] is not None]
            if errors:
                self.ax.bar(compared_scales[:len(errors)], errors, width=0.5,
                            color=["#fb7185" if e > 0 else "#38bdf8" for e in errors])
                self.ax.axhline(0, color="gray")
        elif graph == GRAPH_EFFECTIVE_T:
            self._style_axis(graph, scale_label, "V/V0")
            self.ax.plot(scales, [r["predicted_transmission"] for r in rows], marker="o", color="#38bdf8", label="Predicted")
            if compared:
                self.ax.scatter(compared_scales, [r["measured_transmission"] for r in compared],
                                color="#fb7185", s=60, zorder=3, label="Measured")
            self._legend()
        else:
            self._style_axis(graph, scale_label, "ln(V0/V)")
            self.ax.plot(scales, [r["predicted_attenuation_natural"] for r in rows], marker="o", color="#a78bfa", label="Predicted")
            if compared:
                self.ax.scatter(compared_scales, [r["measured_attenuation_natural"] for r in compared],
                                color="#fb7185", s=60, zorder=3, label="Measured")
            self._legend()
        self.fig.tight_layout()
        self.canvas.draw_idle()

    def _metadata(self):
        spectral = self.analysis["spectral"]
        material = self.materials[self.material_var.get()]
        detector = self.detectors[self.detector_var.get()]
        row = self._selected_row()
        return {
            "material": self.material_var.get(), "spectrum_source": material["source"],
            "detector": detector["display_name"], "weighting_mode": self.weighting_var.get(),
            "wavelength_range": f"{spectral['integration_range_nm'][0]:g}–{spectral['integration_range_nm'][1]:g}",
            "no_film_voltage_mv": float(self.no_film_var.get()),
            "simulation_input": self._scale_label(),
            "inspected_scale": self.analysis["selected_scale"],
            "inspected_measured_mv": row["measured_voltage_mv"],
            "inspected_predicted_mv": row["predicted_voltage_mv"],
            "inspected_error_mv": row["error_mv"],
            "inspected_error_percent": row["percent_error"],
            "reference_thickness": self.reference_thickness_var.get(),
            "absorbance_mode": self.mode_var.get(),
            "baseline_correction": self.baseline_var.get(),
            "baseline_offset": self.analysis["baseline_offset"],
            "clamp_negative": self.clamp_var.get(),
            "fresnel": f"Yes, n = {self.refractive_index_var.get()}" if self.fresnel_var.get() else "No",
            "layer_interface_transmission": self.analysis["layer_interface_transmission"],
            "notes": material["note"] + " " + detector["responsivity_note"],
        }

    def export_excel(self):
        if self.analysis is None:
            return
        default_name = f"broadband_{self.material_var.get()}_{self.detector_var.get()}.xlsx"
        default_name = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in default_name)
        path = filedialog.asksaveasfilename(parent=self, title="Export broadband analysis", initialfile=default_name,
                                            defaultextension=".xlsx", filetypes=[("Excel workbook", "*.xlsx")])
        if not path:
            return
        try:
            experimental_rows = []
            for detector in self.detectors.values():
                for material, measurements in detector.get("materials", {}).items():
                    for layer, voltage in measurements.items():
                        experimental_rows.append({"material": material, "detector": detector["display_name"],
                                                  "layer_count": int(layer), "measured_voltage_mv": voltage,
                                                  "no_film_voltage_mv": detector["no_film_voltage_mv"]})
            saved = export_broadband_workbook(path, self.analysis, self._metadata(), experimental_rows)
            messagebox.showinfo("Export complete", f"Workbook saved to:\n{saved}", parent=self)
        except Exception as exc:
            messagebox.showerror("Export failed", str(exc), parent=self)
