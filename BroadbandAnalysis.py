"""Broadband Beer-Lambert physics primitives.

This module has no UI, file-loading or plotting dependencies.  Wavelengths are
nanometres and voltages are millivolts.  A weighting curve (source spectrum,
detector responsivity, filter/window transmission) may be supplied as a
``(wavelength_nm, values)`` pair or a mapping with ``wavelength_nm`` and
``values`` keys.  Every weighting curve is treated as ZERO outside its sampled
wavelength range; nothing is extrapolated.

The forward model for one film is::

    T(lambda, x)   = 10^(-A_ref(lambda) * x / x_ref)          (base-10)
                   = exp(-A_ref(lambda) * x / x_ref)          (natural)
    W(lambda)      = S(lambda) * R(lambda) * prod(F_i(lambda))
    T_spectral     = int W T dlambda / int W dlambda           (covered range)
    T_effective    = T_spectral * T_interface
    V_pred         = V_dark + (V0 - V_dark) * T_effective
"""

from __future__ import annotations

import math
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np

# Below this fraction of the detector/source weight lying inside the material
# spectrum, a prediction is flagged as PARTIAL SPECTRAL COVERAGE.
DEFAULT_COVERAGE_THRESHOLD = 0.99


def _integrate(values, wavelength_nm):
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is None:  # NumPy < 2.0 compatibility.
        trapezoid = np.trapz
    values = np.asarray(values, dtype=float)
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    if values.size < 2:
        return 0.0
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
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    unique_x, unique_indices = np.unique(x, return_index=True)
    return unique_x, y[unique_indices]


def _curve_xy(curve, name):
    """Parse a weighting curve.  Bare arrays are rejected on purpose: pairing
    values with an implicit grid silently misaligns data that was not already
    sorted and de-duplicated."""
    if curve is None:
        return None
    if isinstance(curve, Mapping):
        curve_x = curve.get("wavelength_nm")
        curve_y = curve.get("values", curve.get("value"))
        if curve_x is None or curve_y is None:
            raise ValueError(f"{name} mapping needs wavelength_nm and values.")
        return _clean_xy(curve_x, curve_y)
    if isinstance(curve, (tuple, list)) and len(curve) == 2:
        return _clean_xy(curve[0], curve[1])
    raise ValueError(f"{name} must be a (wavelength_nm, values) pair or a mapping.")


def absorbance_to_transmission(absorbance, thickness_scale=1.0, absorbance_mode="base10",
                               clamp_negative=False):
    """Convert reference absorbance to transmission at ``thickness_scale = x / x_ref``.

    ``base10``: A = -log10(T)  ->  T = 10^(-A * scale)
    ``natural``: A = -ln(T)    ->  T = exp(-A * scale)
    """
    a = np.asarray(absorbance, dtype=float)
    scale = float(thickness_scale)
    if scale < 0:
        raise ValueError("Thickness scale x / x_ref must be non-negative.")
    if clamp_negative:
        a = np.maximum(a, 0.0)
    mode = str(absorbance_mode).strip().lower().replace("-", "")
    if mode in {"base10", "decadic", "log10"}:
        return np.power(10.0, -(a * scale))
    if mode in {"natural", "napierian", "ln"}:
        return np.exp(-(a * scale))
    raise ValueError("absorbance_mode must be 'base10' or 'natural'.")


def thickness_scale(thickness_um, reference_thickness_um):
    """x / x_ref, requiring both physical thicknesses."""
    if reference_thickness_um is None or not float(reference_thickness_um) > 0:
        raise ValueError("A positive reference-sample thickness x_ref (µm) is required for physical thickness.")
    if thickness_um is None or float(thickness_um) < 0:
        raise ValueError("Thickness must be a non-negative number of µm.")
    return float(thickness_um) / float(reference_thickness_um)


def estimate_baseline_offset(wavelength_nm, absorbance, wavelength_min=None, wavelength_max=None,
                             percentile=5.0):
    """Estimate a flat absorbance baseline offset inside the detector band.

    Heuristic: assumes at least ``percentile`` % of the in-band points do not
    absorb.  Optional correction only; never applied by default.
    """
    x, a = _clean_xy(wavelength_nm, absorbance)
    band = np.ones_like(x, dtype=bool)
    if wavelength_min is not None and np.isfinite(wavelength_min):
        band &= x >= float(wavelength_min)
    if wavelength_max is not None and np.isfinite(wavelength_max):
        band &= x <= float(wavelength_max)
    if not np.any(band):
        raise ValueError("No absorbance points inside the baseline estimation band.")
    return float(np.percentile(a[band], percentile))


def fresnel_layer_transmission(refractive_index, surrounding_index=1.0):
    """Single-pass transmission through the two surrounding/film interfaces of one layer.

    Normal incidence, non-absorbing interfaces, multiple reflections ignored:
    T = (1 - R)^2 with R = ((n1 - n2) / (n1 + n2))^2.
    """
    n1, n2 = float(surrounding_index), float(refractive_index)
    if n1 <= 0 or n2 <= 0:
        raise ValueError("Refractive indices must be positive.")
    reflectance = ((n1 - n2) / (n1 + n2)) ** 2
    return (1.0 - reflectance) ** 2


def blackbody_relative_spectrum(wavelength_nm, temperature_k):
    """Planck spectral radiance per unit wavelength, normalized to its maximum on the grid."""
    temperature_k = float(temperature_k)
    if temperature_k <= 0:
        raise ValueError("Blackbody temperature must be positive.")
    h, c, k = 6.62607015e-34, 299792458.0, 1.380649e-23
    wl_m = np.asarray(wavelength_nm, dtype=float) * 1e-9
    exponent = np.clip(h * c / (wl_m * k * temperature_k), 1e-12, 700)
    intensity = (2 * h * c ** 2 / wl_m ** 5) / np.expm1(exponent)
    maximum = np.nanmax(intensity)
    return intensity / maximum if maximum > 0 else np.zeros_like(wl_m)


def combine_material_absorbances(components, absorbance_mode="base10", clamp_negative=False):
    """Total transmission of a stack of different materials on one shared grid.

    Each component is a mapping with ``absorbance`` (reference spectrum on the
    shared grid) and ``thickness_scale`` (x_i / x_ref_i):
    T_total = 10^(-sum A_i * x_i / x_ref_i).
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


def weighting_domain(source=None, responsivity=None, multiplicative_terms=None,
                     wavelength_min=None, wavelength_max=None, fallback_range=None):
    """Wavelength interval where W(lambda) can be non-zero.

    It is the explicit window intersected with the sampled support of every
    weighting curve (each curve is zero outside its support).
    """
    curves = [_curve_xy(curve, name) for name, curve in
              [("source", source), ("responsivity", responsivity)]
              + [(f"term_{i}", c) for i, c in enumerate(multiplicative_terms or [])]]
    curves = [curve for curve in curves if curve is not None]
    lower = -math.inf if wavelength_min is None else float(wavelength_min)
    upper = math.inf if wavelength_max is None else float(wavelength_max)
    for curve_x, _curve_y in curves:
        lower, upper = max(lower, curve_x[0]), min(upper, curve_x[-1])
    if not np.isfinite(lower) or not np.isfinite(upper):
        if fallback_range is None:
            raise ValueError("Flat weighting needs an explicit wavelength window.")
        lower = fallback_range[0] if not np.isfinite(lower) else lower
        upper = fallback_range[1] if not np.isfinite(upper) else upper
    if upper <= lower:
        raise ValueError("The weighting curves and wavelength window do not overlap.")
    return float(lower), float(upper)


def calculate_weighted_transmission(wavelength_nm, transmission, source=None, responsivity=None,
                                    multiplicative_terms: Optional[Iterable] = None,
                                    wavelength_min=None, wavelength_max=None,
                                    coverage_threshold=DEFAULT_COVERAGE_THRESHOLD):
    """Integrate W(lambda) T(lambda) over the common wavelength grid.

    * The weight domain is the window intersected with every curve's support.
    * Every curve is interpolated onto one common grid; outside its sampled
      range a curve is zero.
    * T is only known where the material spectrum exists.  The effective
      transmission is normalized over that covered part of W, and
      ``coverage_fraction`` reports how much of W that is.  Bounds assume the
      uncovered weight transmits nothing (lower) or everything (upper).
    """
    material_x, material_t = _clean_xy(wavelength_nm, transmission)
    source_xy = _curve_xy(source, "source")
    response_xy = _curve_xy(responsivity, "responsivity")
    term_xy = [_curve_xy(curve, f"term_{i}") for i, curve in enumerate(multiplicative_terms or [])]
    lower, upper = weighting_domain(source_xy, response_xy, term_xy, wavelength_min, wavelength_max,
                                    fallback_range=(material_x[0], material_x[-1]))
    covered_lower, covered_upper = max(lower, material_x[0]), min(upper, material_x[-1])
    if covered_upper <= covered_lower:
        raise ValueError("The material spectrum does not overlap the detector/source weighting range.")

    pieces = [np.array([lower, upper, covered_lower, covered_upper]),
              material_x[(material_x >= lower) & (material_x <= upper)]]
    for curve in [source_xy, response_xy, *term_xy]:
        if curve is not None:
            pieces.append(curve[0][(curve[0] >= lower) & (curve[0] <= upper)])
    grid = np.unique(np.concatenate(pieces))

    def on_grid(curve):
        return None if curve is None else np.interp(grid, curve[0], curve[1], left=0.0, right=0.0)

    mapped_source, mapped_response = on_grid(source_xy), on_grid(response_xy)
    mapped_terms = [on_grid(curve) for curve in term_xy]
    weight = np.ones_like(grid)
    for values in [mapped_source, mapped_response, *mapped_terms]:
        if values is not None:
            weight = weight * values

    covered = (grid >= covered_lower) & (grid <= covered_upper)
    mapped_t = np.full_like(grid, np.nan)
    mapped_t[covered] = np.interp(grid[covered], material_x, material_t)
    weighted_transmitted = weight * mapped_t

    weight_total = _integrate(weight, grid)
    weight_covered = _integrate(weight[covered], grid[covered])
    if not np.isfinite(weight_total) or weight_total <= 0:
        raise ValueError("The integrated detector/source weight must be greater than zero.")
    if weight_covered <= 0:
        raise ValueError("No detector/source weight lies inside the material spectrum range.")
    transmitted_covered = _integrate(weighted_transmitted[covered], grid[covered])
    coverage = weight_covered / weight_total
    uncovered = weight_total - weight_covered
    return {
        "wavelength_nm": grid,
        "transmission": mapped_t,
        "source": mapped_source,
        "responsivity": mapped_response,
        "multiplicative_terms": mapped_terms,
        "weight": weight,
        "weighted_transmitted": weighted_transmitted,
        "covered_mask": covered,
        "weight_integral_total": weight_total,
        "weight_integral_covered": weight_covered,
        "weighted_transmitted_integral": transmitted_covered,
        "coverage_fraction": float(coverage),
        "partial_coverage": bool(coverage < coverage_threshold),
        "effective_transmission": float(transmitted_covered / weight_covered),
        "effective_transmission_bounds": (float(transmitted_covered / weight_total),
                                          float((transmitted_covered + uncovered) / weight_total)),
        "integration_range_nm": (float(lower), float(upper)),
        "material_coverage_nm": (float(covered_lower), float(covered_upper)),
    }


def dark_corrected_voltage(transmission, no_film_voltage_mv, dark_voltage_mv=0.0):
    """V = V_dark + (V0 - V_dark) * T."""
    v0, vd = float(no_film_voltage_mv), float(dark_voltage_mv)
    if v0 <= vd:
        raise ValueError("No-film voltage V0 must be greater than dark voltage V_dark.")
    voltage = vd + (v0 - vd) * np.asarray(transmission, dtype=float)
    return float(voltage) if voltage.ndim == 0 else voltage


def measured_transmission(voltage_mv, no_film_voltage_mv, dark_voltage_mv=0.0):
    """T = (V - V_dark) / (V0 - V_dark)."""
    v0, vd = float(no_film_voltage_mv), float(dark_voltage_mv)
    if v0 <= vd:
        raise ValueError("No-film voltage V0 must be greater than dark voltage V_dark.")
    return (float(voltage_mv) - vd) / (v0 - vd)


def predict_detector_voltage(wavelength_nm, absorbance, no_film_voltage_mv, thickness_scale=1.0,
                             source=None, responsivity=None, multiplicative_terms=None,
                             wavelength_min=None, wavelength_max=None, absorbance_mode="base10",
                             clamp_negative=False, interface_transmission=1.0, dark_voltage_mv=0.0,
                             coverage_threshold=DEFAULT_COVERAGE_THRESHOLD):
    """Predict detector voltage for one thickness scale x / x_ref.

    ``interface_transmission`` is the total wavelength-independent interface
    factor for the whole stack (already raised to the number of layers).
    """
    x, a = _clean_xy(wavelength_nm, absorbance)
    used_absorbance = np.maximum(a, 0.0) if clamp_negative else a
    transmission = absorbance_to_transmission(used_absorbance, thickness_scale, absorbance_mode)
    result = calculate_weighted_transmission(x, transmission, source, responsivity, multiplicative_terms,
                                             wavelength_min, wavelength_max, coverage_threshold)
    covered = result["covered_mask"]
    absorbance_on_grid = np.full_like(result["wavelength_nm"], np.nan)
    absorbance_on_grid[covered] = np.interp(result["wavelength_nm"][covered], x, used_absorbance)
    interface = float(interface_transmission)
    spectral = result["effective_transmission"]
    low, high = result["effective_transmission_bounds"]
    result.update({
        "absorbance": absorbance_on_grid,
        "thickness_scale": float(thickness_scale),
        "spectral_transmission": spectral,
        "interface_transmission": interface,
        "effective_transmission": spectral * interface,
        "effective_transmission_bounds": (low * interface, high * interface),
        "dark_voltage_mv": float(dark_voltage_mv),
        "no_film_voltage_mv": float(no_film_voltage_mv),
        "predicted_voltage_mv": dark_corrected_voltage(spectral * interface, no_film_voltage_mv, dark_voltage_mv),
        "predicted_voltage_bounds_mv": (dark_corrected_voltage(low * interface, no_film_voltage_mv, dark_voltage_mv),
                                        dark_corrected_voltage(high * interface, no_film_voltage_mv, dark_voltage_mv)),
    })
    return result


def agreement_metrics(measured: Sequence[float], predicted: Sequence[float]):
    """MAE, RMSE, MAPE and agreement R^2 = 1 - SS_res/SS_tot against the 1:1 line."""
    measured = np.asarray(measured, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if measured.size == 0:
        return {"n": 0, "mae_mv": None, "rmse_mv": None, "mape_percent": None, "r_squared": None}
    residual = predicted - measured
    nonzero = measured != 0
    metrics = {
        "n": int(measured.size),
        "mae_mv": float(np.mean(np.abs(residual))),
        "rmse_mv": float(np.sqrt(np.mean(residual ** 2))),
        "mape_percent": float(np.mean(np.abs(residual[nonzero] / measured[nonzero])) * 100) if np.any(nonzero) else None,
        "r_squared": None,
    }
    denominator = float(np.sum((measured - np.mean(measured)) ** 2))
    if measured.size >= 2 and denominator > 0:
        metrics["r_squared"] = float(1.0 - np.sum(residual ** 2) / denominator)
    return metrics


# Bisection is used instead of a closed-form inverse because T_effective is a
# weighted integral of 10^(-A(lambda)*scale): it has no analytic inverse, but it
# IS strictly decreasing in scale, so a bracket plus bisection is exact to
# machine precision and cannot converge to a spurious root.
DEFAULT_MAX_THICKNESS_SCALE = 1000.0


def solve_thickness_scale(transmission_of_scale, target_transmission,
                          max_scale=DEFAULT_MAX_THICKNESS_SCALE, tolerance=1e-10, max_iterations=200):
    """Invert a monotonically decreasing T(scale) for the thickness scale x / x_ref.

    ``transmission_of_scale`` must be the SAME forward model used for
    prediction, evaluated at an arbitrary (fractional) scale.  Returns a
    mapping with ``thickness_scale``, the transmission actually reached, the
    final bracket and a ``status``:

    * ``ok`` - a bracket was found and bisected to ``tolerance``.
    * ``at_or_above_no_film`` - the target transmission is >= T(0); the film
      cannot be thinner than zero, so the scale is reported as 0.
    * ``exceeds_max_scale`` - even at ``max_scale`` the model transmits more
      than the target; the scale is a lower bound, not an estimate.
    """
    target = float(target_transmission)
    if not math.isfinite(target):
        raise ValueError("Target transmission must be a finite number.")
    if target <= 0:
        raise ValueError("Target transmission must be greater than zero "
                         "(the measured voltage is at or below the dark voltage).")
    if max_scale <= 0:
        raise ValueError("The maximum thickness scale must be positive.")

    def evaluate(scale):
        # Absorbance that is negative anywhere makes T(lambda) = 10^(-A*scale) overflow at
        # large scales; an overflow means "transmits far more than the target", never a root.
        value = float(transmission_of_scale(scale))
        return math.inf if not math.isfinite(value) else value

    t_zero = evaluate(0.0)
    if target >= t_zero:
        return {"thickness_scale": 0.0, "transmission": t_zero, "target_transmission": target,
                "status": "at_or_above_no_film", "converged": False, "iterations": 0, "bracket": (0.0, 0.0)}

    lower, upper = 0.0, 1.0
    t_upper = evaluate(upper)
    while t_upper > target:
        if upper >= max_scale:
            return {"thickness_scale": float(max_scale), "transmission": t_upper, "target_transmission": target,
                    "status": "exceeds_max_scale", "converged": False, "iterations": 0,
                    "bracket": (float(max_scale), math.inf)}
        lower = upper
        upper = min(upper * 2.0, max_scale)
        t_upper = evaluate(upper)

    iterations = 0
    while iterations < max_iterations and (upper - lower) > tolerance * max(1.0, upper):
        middle = 0.5 * (lower + upper)
        if evaluate(middle) > target:
            lower = middle
        else:
            upper = middle
        iterations += 1
    scale = 0.5 * (lower + upper)
    return {"thickness_scale": float(scale), "transmission": evaluate(scale),
            "target_transmission": target, "status": "ok",
            "converged": bool((upper - lower) <= tolerance * max(1.0, upper)),
            "iterations": iterations, "bracket": (float(lower), float(upper))}
