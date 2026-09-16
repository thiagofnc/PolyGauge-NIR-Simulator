"""Broadband Beer-Lambert calculations for absorbance spectra.

This module deliberately has no UI dependencies.  Wavelengths are nanometres,
voltages are millivolts, and a spectral curve may be supplied as either an
array on the material grid, a ``(wavelength_nm, values)`` pair, or a mapping
with ``wavelength_nm`` and ``values`` keys.
"""

from __future__ import annotations

import math
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np


def _integrate(values, wavelength_nm):
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is None:  # NumPy < 2.0 compatibility.
        trapezoid = np.trapz
    return float(trapezoid(values, wavelength_nm))


def _clean_xy(wavelength_nm, values):
    x = np.asarray(wavelength_nm, dtype=float).reshape(-1)
    y = np.asarray(values, dtype=float).reshape(-1)
    if x.size != y.size:
        raise ValueError("Wavelength and spectral value arrays must have the same length.")
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if x.size < 2:
        raise ValueError("A spectrum requires at least two finite wavelength points.")
    order = np.argsort(x)
    x, y = x[order], y[order]
    unique_x, unique_indices = np.unique(x, return_index=True)
    return unique_x, y[unique_indices]


def _curve_xy(curve, default_wavelength_nm, name):
    if curve is None:
        return None
    if isinstance(curve, Mapping):
        curve_x = curve.get("wavelength_nm")
        curve_y = curve.get("values", curve.get("value"))
        if curve_x is None or curve_y is None:
            raise ValueError(f"{name} mapping needs wavelength_nm and values.")
        return _clean_xy(curve_x, curve_y)
    if isinstance(curve, (tuple, list)) and len(curve) == 2:
        first = np.asarray(curve[0])
        second = np.asarray(curve[1])
        if first.ndim > 0 and second.ndim > 0:
            return _clean_xy(first, second)
    values = np.asarray(curve, dtype=float)
    return _clean_xy(default_wavelength_nm, values)


def absorbance_to_transmission(
    absorbance,
    thickness_scale=1.0,
    absorbance_mode="base10",
    clamp_negative=False,
):
    """Convert absorbance/attenuation to transmission.

    ``base10`` means spectroscopic absorbance A=-log10(T). ``natural`` and
    ``napierian`` mean A=-ln(T).  ``alpha`` is accepted as a synonym when the
    caller has already folded physical thickness into ``thickness_scale``.
    """
    a = np.asarray(absorbance, dtype=float)
    scale = np.asarray(thickness_scale, dtype=float)
    if clamp_negative:
        a = np.maximum(a, 0.0)
    mode = str(absorbance_mode).strip().lower().replace("-", "")
    if mode in {"base10", "decadic", "log10"}:
        return np.power(10.0, -(a * scale))
    if mode in {"natural", "napierian", "ln", "alpha"}:
        return np.exp(-(a * scale))
    raise ValueError("absorbance_mode must be 'base10' or 'natural'.")


def estimate_baseline_offset(wavelength_nm, absorbance, wavelength_min=None, wavelength_max=None,
                             percentile=5.0):
    """Estimate a flat absorbance baseline offset inside the detector band.

    FTIR scans whose background does not exactly match the sample often sit at
    a constant non-zero absorbance (e.g. -0.03 A, i.e. ~107 %T).  Integrated over
    a broad band that offset swamps the real absorption peaks, so a low
    percentile of the in-band absorbance is treated as the zero line.
    """
    x, a = _clean_xy(wavelength_nm, absorbance)
    band = np.ones_like(x, dtype=bool)
    if wavelength_min is not None:
        band &= x >= float(wavelength_min)
    if wavelength_max is not None:
        band &= x <= float(wavelength_max)
    if not np.any(band):
        band = np.ones_like(x, dtype=bool)
    return float(np.percentile(a[band], percentile))


def fresnel_layer_transmission(refractive_index, surrounding_index=1.0):
    """Single-pass transmission through the two air/film interfaces of one layer.

    Normal incidence, multiple reflections ignored: T = (1 - R)^2 with
    R = ((n1 - n2) / (n1 + n2))^2.
    """
    n1, n2 = float(surrounding_index), float(refractive_index)
    if n1 <= 0 or n2 <= 0:
        raise ValueError("Refractive indices must be positive.")
    reflectance = ((n1 - n2) / (n1 + n2)) ** 2
    return (1.0 - reflectance) ** 2


def combine_material_absorbances(components, absorbance_mode="base10", clamp_negative=False):
    """Return total transmission for multiple material contributions.

    Each component is a mapping containing ``absorbance`` and optionally
    ``thickness_scale``.  All arrays must already share a wavelength grid.
    This forward model intentionally makes no claim that component thicknesses
    can be recovered from one broadband measurement.
    """
    total_exponent = None
    for component in components:
        values = np.asarray(component["absorbance"], dtype=float)
        if clamp_negative:
            values = np.maximum(values, 0.0)
        contribution = values * float(component.get("thickness_scale", 1.0))
        total_exponent = contribution.copy() if total_exponent is None else total_exponent + contribution
    if total_exponent is None:
        raise ValueError("At least one material component is required.")
    return absorbance_to_transmission(total_exponent, absorbance_mode=absorbance_mode)


def calculate_weighted_transmission(
    wavelength_nm,
    transmission,
    source=None,
    responsivity=None,
    multiplicative_terms: Optional[Iterable] = None,
    wavelength_min=None,
    wavelength_max=None,
):
    """Integrate W(lambda)T(lambda)/W(lambda) on the valid common overlap.

    The common grid contains every input sample inside the overlap. Detector
    response is never extrapolated: providing a response curve restricts the
    integration interval to that curve's measured/support range.
    """
    material_x, material_t = _clean_xy(wavelength_nm, transmission)
    named_curves = [("source", source), ("responsivity", responsivity)]
    named_curves.extend((f"term_{i}", curve) for i, curve in enumerate(multiplicative_terms or []))
    parsed = [(name, _curve_xy(curve, material_x, name)) for name, curve in named_curves]
    parsed = [(name, curve) for name, curve in parsed if curve is not None]

    lower = material_x[0] if wavelength_min is None else max(material_x[0], float(wavelength_min))
    upper = material_x[-1] if wavelength_max is None else min(material_x[-1], float(wavelength_max))
    for _name, (curve_x, _curve_y) in parsed:
        lower, upper = max(lower, curve_x[0]), min(upper, curve_x[-1])
    if not np.isfinite(lower) or not np.isfinite(upper) or upper <= lower:
        raise ValueError("The required spectra and wavelength range do not overlap.")

    pieces = [material_x[(material_x >= lower) & (material_x <= upper)], np.array([lower, upper])]
    for _name, (curve_x, _curve_y) in parsed:
        pieces.append(curve_x[(curve_x >= lower) & (curve_x <= upper)])
    common_x = np.unique(np.concatenate(pieces))
    mapped_t = np.interp(common_x, material_x, material_t)
    mapped = {name: np.interp(common_x, x, y) for name, (x, y) in parsed}

    weight = np.ones_like(common_x)
    for values in mapped.values():
        weight *= values
    denominator = _integrate(weight, common_x)
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError("The integrated spectral weight must be greater than zero.")
    weighted_transmitted = weight * mapped_t
    effective = _integrate(weighted_transmitted, common_x) / denominator
    return {
        "wavelength_nm": common_x,
        "transmission": mapped_t,
        "source": mapped.get("source"),
        "responsivity": mapped.get("responsivity"),
        "multiplicative_terms": {k: v for k, v in mapped.items() if k.startswith("term_")},
        "weight": weight,
        "weighted_input": weight.copy(),
        "weighted_transmitted": weighted_transmitted,
        "effective_transmission": float(effective),
        "integration_range_nm": (float(lower), float(upper)),
    }


def predict_detector_voltage(
    wavelength_nm,
    absorbance,
    no_film_voltage_mv,
    thickness_scale=1.0,
    source=None,
    responsivity=None,
    multiplicative_terms=None,
    wavelength_min=None,
    wavelength_max=None,
    absorbance_mode="base10",
    clamp_negative=False,
    layer_interface_transmission=1.0,
):
    """Predict a normalized broadband detector voltage and return spectra.

    ``layer_interface_transmission`` is a wavelength-independent loss applied
    once per layer (e.g. :func:`fresnel_layer_transmission`), so the reported
    ``effective_transmission`` is ``spectral_transmission * interface**scale``.
    """
    x, a = _clean_xy(wavelength_nm, absorbance)
    transmission = absorbance_to_transmission(a, thickness_scale, absorbance_mode, clamp_negative)
    result = calculate_weighted_transmission(
        x, transmission, source, responsivity, multiplicative_terms,
        wavelength_min, wavelength_max,
    )
    interface = float(layer_interface_transmission) ** float(thickness_scale)
    used_absorbance = np.maximum(a, 0.0) if clamp_negative else a
    result["absorbance"] = np.interp(result["wavelength_nm"], x, used_absorbance)
    result["thickness_scale"] = float(thickness_scale)
    result["spectral_transmission"] = result["effective_transmission"]
    result["interface_transmission"] = interface
    result["effective_transmission"] = result["spectral_transmission"] * interface
    result["predicted_voltage_mv"] = float(no_film_voltage_mv) * result["effective_transmission"]
    return result


def _agreement_metrics(measured, predicted):
    measured = np.asarray(measured, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    residual = predicted - measured
    metrics = {
        "mae_mv": float(np.mean(np.abs(residual))),
        "rmse_mv": float(np.sqrt(np.mean(residual ** 2))),
        "r_squared": None,
    }
    denominator = float(np.sum((measured - np.mean(measured)) ** 2))
    if measured.size >= 2 and denominator > 0:
        metrics["r_squared"] = float(1.0 - np.sum(residual ** 2) / denominator)
    return metrics


def analyze_layers(
    wavelength_nm,
    absorbance,
    no_film_voltage_mv,
    thickness_scales: Sequence[float],
    measured_voltage_by_scale: Optional[Mapping[float, float]] = None,
    selected_scale=None,
    **prediction_options,
):
    """Calculate every requested layer/thickness scale from one shared model."""
    scales = [float(value) for value in thickness_scales]
    if not scales:
        raise ValueError("At least one thickness scale is required.")
    measured_voltage_by_scale = measured_voltage_by_scale or {}
    rows, spectra_by_scale = [], {}
    for scale in scales:
        spectral = predict_detector_voltage(
            wavelength_nm, absorbance, no_film_voltage_mv, scale, **prediction_options
        )
        spectra_by_scale[scale] = spectral
        predicted_v = spectral["predicted_voltage_mv"]
        predicted_t = spectral["effective_transmission"]
        measured_v = measured_voltage_by_scale.get(scale)
        if measured_v is None and scale.is_integer():
            measured_v = measured_voltage_by_scale.get(int(scale))
        row = {
            "layer_count": int(scale) if scale.is_integer() else scale,
            "thickness_scale": scale,
            "measured_voltage_mv": None if measured_v is None else float(measured_v),
            "predicted_voltage_mv": predicted_v,
            "measured_transmission": None,
            "predicted_transmission": predicted_t,
            "spectral_transmission": spectral["spectral_transmission"],
            "interface_transmission": spectral["interface_transmission"],
            "measured_attenuation_natural": None,
            "predicted_attenuation_natural": -math.log(predicted_t) if predicted_t > 0 else math.inf,
            "predicted_absorbance_base10": -math.log10(predicted_t) if predicted_t > 0 else math.inf,
            "error_mv": None,
            "absolute_error_mv": None,
            "percent_error": None,
            "absolute_percent_error": None,
        }
        if measured_v is not None:
            measured_t = float(measured_v) / float(no_film_voltage_mv)
            error = predicted_v - float(measured_v)
            percent = 100.0 * error / float(measured_v) if measured_v != 0 else None
            row.update({
                "measured_transmission": measured_t,
                "measured_attenuation_natural": -math.log(measured_t) if measured_t > 0 else math.inf,
                "error_mv": error,
                "absolute_error_mv": abs(error),
                "percent_error": percent,
                "absolute_percent_error": None if percent is None else abs(percent),
            })
        rows.append(row)

    compared = [row for row in rows if row["measured_voltage_mv"] is not None]
    metrics = _agreement_metrics(
        [row["measured_voltage_mv"] for row in compared],
        [row["predicted_voltage_mv"] for row in compared],
    ) if compared else {"mae_mv": None, "rmse_mv": None, "r_squared": None}
    percent_errors = [row["absolute_percent_error"] for row in compared
                      if row["absolute_percent_error"] is not None]
    metrics["mape_percent"] = float(np.mean(percent_errors)) if percent_errors else None
    selected = float(selected_scale) if selected_scale is not None else scales[0]
    selected = min(scales, key=lambda value: abs(value - selected))
    return {
        "rows": rows,
        "metrics": metrics,
        "spectral": spectra_by_scale[selected],
        "spectra_by_scale": spectra_by_scale,
        "selected_scale": selected,
        "negative_absorbance_present": bool(np.any(np.asarray(absorbance, dtype=float) < 0)),
    }
