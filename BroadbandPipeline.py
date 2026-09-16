"""Single shared broadband analysis pipeline.

``run_broadband_analysis(config)`` is the only place where a material spectrum
becomes a predicted detector voltage.  The UI, plots, Excel export and any
future thickness solver consume its result and must not recompute physics.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from BroadbandAnalysis import (DEFAULT_COVERAGE_THRESHOLD, agreement_metrics, blackbody_relative_spectrum,
                               estimate_baseline_offset, fresnel_layer_transmission, measured_transmission,
                               predict_detector_voltage, weighting_domain)

# Thickness input modes
LAYERS_PHYSICAL = "layers_physical"          # x = n * t_layer, scale = x / x_ref  (reports µm)
LAYERS_RATIO = "layers_ratio"                # scale = n * (t_layer / x_ref), x_ref unknown (no µm)
REFERENCE_MULTIPLIER = "reference_multiplier"  # scale given directly (no µm, no layer matching)

DETECTOR_FLAT_BAND = "flat_band"
DETECTOR_RESPONSIVITY = "responsivity_curve"

SOURCE_FLAT = "flat"
SOURCE_BLACKBODY = "blackbody"
SOURCE_MEASURED = "measured"

BLACKBODY_GRID_POINTS = 4001
MATCH_TOLERANCE = 1e-9


@dataclass
class BroadbandConfig:
    material: str
    values: Sequence[float]
    no_film_voltage_mv: float
    thickness_mode: str = LAYERS_RATIO
    reference_thickness_um: Optional[float] = None      # overrides material metadata when given
    layer_thickness_um: Optional[float] = None          # LAYERS_PHYSICAL
    layer_to_reference_ratio: Optional[float] = None    # LAYERS_RATIO
    dark_voltage_mv: float = 0.0
    wavelength_min_nm: Optional[float] = None
    wavelength_max_nm: Optional[float] = None
    detector_name: str = ""
    detector_mode: str = DETECTOR_FLAT_BAND
    detector_responsivity: Optional[dict] = None        # {"wavelength_nm", "values", "note"}
    source_mode: str = SOURCE_FLAT
    source_label: str = ""
    source_temperature_k: Optional[float] = None
    source_min_nm: Optional[float] = None
    source_max_nm: Optional[float] = None
    source_spectrum: Optional[dict] = None              # {"wavelength_nm", "values", "note"}
    optical_terms: List[dict] = field(default_factory=list)
    absorbance_convention: Optional[str] = None          # None = material metadata
    baseline_correction: bool = False
    clamp_negative: bool = False
    interface_correction: bool = False
    refractive_index: Optional[float] = None
    measured_voltages_by_layer: Optional[Dict[int, float]] = None
    inspect_value: Optional[float] = None
    coverage_threshold: float = DEFAULT_COVERAGE_THRESHOLD


def _curve(spec, name):
    if spec is None:
        raise ValueError(f"{name} data is required for the selected mode.")
    return (np.asarray(spec["wavelength_nm"], dtype=float), np.asarray(spec["values"], dtype=float))


def _weighting(config, material):
    """Build (source, responsivity, terms, window, descriptions) without physics."""
    window_min, window_max = config.wavelength_min_nm, config.wavelength_max_nm
    if config.detector_mode == DETECTOR_FLAT_BAND:
        if window_min is None or window_max is None:
            raise ValueError("Flat band approximation needs both minimum and maximum wavelength.")
        responsivity = None
        detector_desc = f"Flat band approximation R(λ)=1 over {window_min:g}-{window_max:g} nm (NOT a measured responsivity)"
    elif config.detector_mode == DETECTOR_RESPONSIVITY:
        responsivity = _curve(config.detector_responsivity, "Detector responsivity")
        detector_desc = "Loaded detector responsivity: " + config.detector_responsivity.get("note", "")
    else:
        raise ValueError(f"Unknown detector mode {config.detector_mode!r}.")

    terms = [_curve(term, "Optical term") for term in config.optical_terms]
    if config.source_mode == SOURCE_FLAT:
        source = None
        source_desc = "Flat source (no source spectral weighting)"
    elif config.source_mode == SOURCE_MEASURED:
        source = _curve(config.source_spectrum, "Measured source spectrum")
        source_desc = "Measured source spectrum: " + config.source_spectrum.get("note", "")
    elif config.source_mode == SOURCE_BLACKBODY:
        if not config.source_temperature_k:
            raise ValueError("Blackbody approximation needs a source temperature.")
        lower, upper = weighting_domain(None, responsivity, terms,
                                        _max_optional(window_min, config.source_min_nm),
                                        _min_optional(window_max, config.source_max_nm),
                                        fallback_range=(material["wavelength_nm"][0], material["wavelength_nm"][-1]))
        grid = np.linspace(lower, upper, BLACKBODY_GRID_POINTS)
        source = (grid, blackbody_relative_spectrum(grid, config.source_temperature_k))
        source_desc = (f"Blackbody approximation {config.source_temperature_k:g} K over {lower:g}-{upper:g} nm"
                       + (f" ({config.source_label})" if config.source_label else ""))
    else:
        raise ValueError(f"Unknown source mode {config.source_mode!r}.")
    return source, responsivity, terms, (window_min, window_max), detector_desc, source_desc


def _max_optional(a, b):
    values = [v for v in (a, b) if v is not None]
    return max(values) if values else None


def _min_optional(a, b):
    values = [v for v in (a, b) if v is not None]
    return min(values) if values else None


def _thickness_rows(config, reference_thickness_um):
    """Return [(layer_count|None, thickness_um|None, scale)] and a description."""
    values = [float(v) for v in config.values]
    if not values:
        raise ValueError("At least one thickness value is required.")
    if any(v < 0 for v in values):
        raise ValueError("Thickness values must be non-negative.")
    if config.thickness_mode == LAYERS_PHYSICAL:
        if reference_thickness_um is None or reference_thickness_um <= 0:
            raise ValueError("Physical thickness mode needs the reference-sample thickness x_ref (µm). "
                             "Use 'Layers × thickness ratio' or 'Reference-sample multiplier' if it is unknown.")
        if not config.layer_thickness_um or config.layer_thickness_um <= 0:
            raise ValueError("Physical thickness mode needs a positive layer thickness (µm).")
        rows = [(v, v * config.layer_thickness_um, v * config.layer_thickness_um / reference_thickness_um) for v in values]
        desc = (f"x = layers × {config.layer_thickness_um:g} µm; scale = x / x_ref with x_ref = {reference_thickness_um:g} µm")
    elif config.thickness_mode == LAYERS_RATIO:
        ratio = config.layer_to_reference_ratio
        if ratio is None or ratio <= 0:
            raise ValueError("Layer thickness ratio t_layer / x_ref must be positive.")
        rows = [(v, None, v * ratio) for v in values]
        desc = (f"scale = layers × {ratio:g} (ASSUMED t_layer / x_ref ratio; x_ref unknown, so no µm are reported)")
    elif config.thickness_mode == REFERENCE_MULTIPLIER:
        rows = [(None, None, v) for v in values]
        desc = "scale = equivalent number of reference samples (no µm, no layer matching)"
    else:
        raise ValueError(f"Unknown thickness mode {config.thickness_mode!r}.")
    for layer_count, _thickness, _scale in rows:
        if layer_count is not None and not float(layer_count).is_integer():
            raise ValueError("Layer counts must be whole numbers; use the multiplier mode for fractional values.")
    return [(None if n is None else int(n), x, s) for n, x, s in rows], desc


def run_broadband_analysis(config: BroadbandConfig, material_library=None, base_dir=None):
    """Run the complete forward model and comparison.  See module docstring."""
    if material_library is None:
        from SpectralData import load_material_library
        material_library = load_material_library(base_dir)
    if config.material not in material_library:
        raise ValueError(f"Unknown material {config.material!r}.")
    material = material_library[config.material]
    warnings: List[str] = []

    # 1-2. Spectrum and absorbance convention.
    wavelength = np.asarray(material["wavelength_nm"], dtype=float)
    raw_absorbance = np.asarray(material["absorbance"], dtype=float)
    convention = config.absorbance_convention or material["convention"]
    if convention not in ("base10", "natural"):
        raise ValueError("Absorbance convention must be 'base10' or 'natural'.")
    if config.absorbance_convention and config.absorbance_convention != material["convention"]:
        warnings.append(f"Absorbance convention overridden to {convention} (material data is {material['convention']}).")

    # 5-8. Weighting curves and domain.
    source, responsivity, terms, window, detector_desc, source_desc = _weighting(config, material)
    domain = weighting_domain(source, responsivity, terms, *window,
                              fallback_range=(wavelength[0], wavelength[-1]))
    in_domain = (wavelength >= domain[0]) & (wavelength <= domain[1])

    # Negative absorbance within the ACTIVE region.
    raw_in_domain = raw_absorbance[in_domain]
    negative = {"points": int(np.sum(raw_in_domain < 0)), "total_points": int(raw_in_domain.size),
                "percent": float(100 * np.mean(raw_in_domain < 0)) if raw_in_domain.size else 0.0,
                "min_absorbance": float(np.min(raw_in_domain)) if raw_in_domain.size else None}
    if negative["points"]:
        warnings.append(f"{negative['points']} of {negative['total_points']} absorbance points "
                        f"({negative['percent']:.1f}%) in the active {domain[0]:.0f}-{domain[1]:.0f} nm region are negative "
                        f"(min {negative['min_absorbance']:.4f}). Uncorrected, these raise transmission above 1.")

    # 3-4. Optional corrections (explicit and reported).
    baseline_offset = 0.0
    if config.baseline_correction:
        baseline_offset = estimate_baseline_offset(wavelength, raw_absorbance, *domain)
    corrected_absorbance = raw_absorbance - baseline_offset
    clamped_points = int(np.sum(corrected_absorbance[in_domain] < 0)) if config.clamp_negative else 0
    interface_per_layer = 1.0
    reference_thickness_um = config.reference_thickness_um or material.get("reference_thickness_um")
    thickness_rows, thickness_desc = _thickness_rows(config, reference_thickness_um)
    if config.interface_correction:
        if config.refractive_index is None:
            raise ValueError("Interface correction needs a refractive index.")
        if config.thickness_mode == REFERENCE_MULTIPLIER:
            raise ValueError("Interface correction is applied per physical layer and needs a layer-count mode.")
        interface_per_layer = fresnel_layer_transmission(config.refractive_index)
    corrections = {
        "baseline_enabled": config.baseline_correction, "baseline_offset": baseline_offset,
        "baseline_method": "5th percentile of in-band absorbance subtracted" if config.baseline_correction else None,
        "clamp_enabled": config.clamp_negative, "clamped_points": clamped_points,
        "clamped_percent": 100.0 * clamped_points / max(1, int(np.sum(in_domain))),
        "interface_enabled": config.interface_correction, "refractive_index": config.refractive_index,
        "interface_per_layer": interface_per_layer,
        "interface_assumptions": ("Normal incidence; two air/film interfaces per layer (air gaps between layers); "
                                  "wavelength-independent n; no multiple reflections or interference."
                                  if config.interface_correction else None),
    }
    corrections["any_enabled"] = bool(config.baseline_correction or config.clamp_negative or config.interface_correction)

    # 9-14. Per thickness: scaling, transmission, weighted integration, interface, dark voltage.
    options = dict(source=source, responsivity=responsivity, multiplicative_terms=terms,
                   wavelength_min=window[0], wavelength_max=window[1], absorbance_mode=convention,
                   dark_voltage_mv=config.dark_voltage_mv, coverage_threshold=config.coverage_threshold)
    measured_by_layer = {int(k): float(v) for k, v in (config.measured_voltages_by_layer or {}).items()}
    if measured_by_layer and config.thickness_mode == REFERENCE_MULTIPLIER:
        warnings.append("Measured data is keyed by layer count and is NOT compared in reference-multiplier mode.")
        measured_by_layer = {}

    rows, raw_spectra, spectra = [], [], []
    for layer_count, thickness_um, scale in thickness_rows:
        raw = predict_detector_voltage(wavelength, raw_absorbance, config.no_film_voltage_mv, scale, **options)
        interface_total = interface_per_layer ** (layer_count or 0)
        corrected = predict_detector_voltage(wavelength, corrected_absorbance, config.no_film_voltage_mv, scale,
                                             clamp_negative=config.clamp_negative,
                                             interface_transmission=interface_total, **options)
        raw_spectra.append(raw)
        spectra.append(corrected)
        row = {
            "layer_count": layer_count, "thickness_um": thickness_um, "thickness_scale": scale,
            "raw_predicted_voltage_mv": raw["predicted_voltage_mv"],
            "raw_effective_transmission": raw["effective_transmission"],
            "predicted_voltage_mv": corrected["predicted_voltage_mv"],
            "predicted_voltage_bounds_mv": corrected["predicted_voltage_bounds_mv"],
            "spectral_transmission": corrected["spectral_transmission"],
            "interface_transmission": corrected["interface_transmission"],
            "predicted_transmission": corrected["effective_transmission"],
            "predicted_attenuation_natural": _neg_log(corrected["effective_transmission"], math.log),
            "predicted_absorbance_base10": _neg_log(corrected["effective_transmission"], math.log10),
            "weight_integral_total": corrected["weight_integral_total"],
            "weight_integral_covered": corrected["weight_integral_covered"],
            "weighted_transmitted_integral": corrected["weighted_transmitted_integral"],
            "coverage_fraction": corrected["coverage_fraction"],
            "partial_coverage": corrected["partial_coverage"],
            "measured_voltage_mv": None, "measured_transmission": None, "measured_attenuation_natural": None,
            "error_mv": None, "absolute_error_mv": None, "percent_error": None, "absolute_percent_error": None,
            "raw_error_mv": None, "raw_percent_error": None,
        }
        # 15. Comparison, only for rows that are genuinely a layer count.
        if layer_count is not None and layer_count in measured_by_layer:
            measured_v = measured_by_layer[layer_count]
            measured_t = measured_transmission(measured_v, config.no_film_voltage_mv, config.dark_voltage_mv)
            error = row["predicted_voltage_mv"] - measured_v
            raw_error = row["raw_predicted_voltage_mv"] - measured_v
            row.update({
                "measured_voltage_mv": measured_v, "measured_transmission": measured_t,
                "measured_attenuation_natural": _neg_log(measured_t, math.log),
                "error_mv": error, "absolute_error_mv": abs(error),
                "percent_error": 100 * error / measured_v if measured_v else None,
                "absolute_percent_error": abs(100 * error / measured_v) if measured_v else None,
                "raw_error_mv": raw_error, "raw_percent_error": 100 * raw_error / measured_v if measured_v else None,
            })
        rows.append(row)

    # 16. Metrics.
    compared = [row for row in rows if row["measured_voltage_mv"] is not None]
    measured_values = [row["measured_voltage_mv"] for row in compared]
    metrics = agreement_metrics(measured_values, [row["predicted_voltage_mv"] for row in compared])
    raw_metrics = agreement_metrics(measured_values, [row["raw_predicted_voltage_mv"] for row in compared])
    coverage_fraction = spectra[0]["coverage_fraction"]
    partial = spectra[0]["partial_coverage"]
    if partial:
        warnings.insert(0, f"PARTIAL SPECTRAL COVERAGE: the material spectrum covers only {coverage_fraction:.1%} of "
                           f"the detector/source weight. Predicted voltages and errors are NOT fully valid; "
                           f"bounds assume the uncovered region transmits 0-100%.")

    inspected_index = None
    if config.inspect_value is not None:
        matches = [i for i, value in enumerate(config.values)
                   if abs(float(value) - float(config.inspect_value)) <= MATCH_TOLERANCE]
        if not matches:
            raise ValueError(f"Inspected value {config.inspect_value:g} is not one of the simulated values.")
        inspected_index = matches[0]

    grid = spectra[0]
    return {
        "config": _config_summary(config),
        "material": {"name": material["name"], "source": material["source"], "note": material["note"],
                     "convention": convention, "reference_thickness_um": reference_thickness_um},
        "thickness": {"mode": config.thickness_mode, "description": thickness_desc,
                      "reference_thickness_um": reference_thickness_um,
                      "layer_thickness_um": config.layer_thickness_um,
                      "layer_to_reference_ratio": config.layer_to_reference_ratio,
                      "reports_um": config.thickness_mode == LAYERS_PHYSICAL},
        "weighting": {"detector_name": config.detector_name, "detector_mode": config.detector_mode,
                      "detector_description": detector_desc, "source_mode": config.source_mode,
                      "source_description": source_desc,
                      "optical_terms": [term.get("note", "") for term in config.optical_terms],
                      "window_nm": window, "weight_domain_nm": domain},
        "voltages": {"no_film_voltage_mv": float(config.no_film_voltage_mv),
                     "dark_voltage_mv": float(config.dark_voltage_mv)},
        "coverage": {"fraction": coverage_fraction, "partial": partial, "threshold": config.coverage_threshold,
                     "material_range_nm": (float(wavelength[0]), float(wavelength[-1])),
                     "covered_range_nm": grid["material_coverage_nm"], "weight_range_nm": grid["integration_range_nm"]},
        "negative_absorbance": negative,
        "corrections": corrections,
        "grid": {"wavelength_nm": grid["wavelength_nm"], "covered_mask": grid["covered_mask"],
                 "raw_absorbance": raw_spectra[0]["absorbance"], "absorbance_used": grid["absorbance"],
                 "source": grid["source"], "responsivity": grid["responsivity"],
                 "optical_terms": grid["multiplicative_terms"],
                 "optical_terms_product": (np.prod(grid["multiplicative_terms"], axis=0)
                                           if grid["multiplicative_terms"] else None),
                 "weight": grid["weight"]},
        "rows": rows, "spectra": spectra, "raw_spectra": raw_spectra,
        "metrics": metrics, "raw_metrics": raw_metrics,
        "inspected_index": inspected_index, "warnings": warnings,
    }


def _neg_log(value, log):
    return -log(value) if value is not None and value > 0 else None


def _config_summary(config):
    summary = asdict(config)
    for key in ("detector_responsivity", "source_spectrum"):
        if summary.get(key):
            summary[key] = summary[key].get("note", "loaded")
    summary["optical_terms"] = [term.get("note", "loaded") for term in config.optical_terms]
    summary["values"] = [float(v) for v in config.values]
    return summary
