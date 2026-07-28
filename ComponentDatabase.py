import json
from pathlib import Path

import numpy as np


def gaussian_peak(wavelengths, center, width, amplitude):
    return amplitude * np.exp(-((wavelengths - center) ** 2) / (2 * width ** 2))


def gaussian_bandpass(wavelengths, center, fwhm, peak=1.0):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return gaussian_peak(wavelengths, center, sigma, peak)


def blackbody_spectrum(wavelengths_nm, temp_k):
    h, c, k = 6.626e-34, 3.0e8, 1.38e-23
    wl_m = wavelengths_nm * 1e-9
    intensity = (2 * h * c**2 / wl_m**5) / (np.exp(h * c / (wl_m * k * temp_k)) - 1)
    peak = np.max(intensity)
    if peak <= 0:
        return np.zeros_like(wavelengths_nm, dtype=float)
    return intensity / peak


def bounded_blackbody_spectrum(wavelengths_nm, temp_k, min_nm=None, max_nm=None):
    spectrum = blackbody_spectrum(wavelengths_nm, temp_k)
    mask = np.ones_like(wavelengths_nm, dtype=bool)
    if min_nm is not None:
        mask &= wavelengths_nm >= min_nm
    if max_nm is not None:
        mask &= wavelengths_nm <= max_nm
    return np.where(mask, spectrum, 0.0)


def piecewise_response(wavelengths_nm, points):
    xs = np.array([point[0] for point in points], dtype=float)
    ys = np.array([point[1] for point in points], dtype=float)
    return np.maximum(np.interp(wavelengths_nm, xs, ys), 0.0)


def default_component_database():
    bundled_database = Path(__file__).with_name("component_database.json")
    if bundled_database.exists():
        with bundled_database.open("r", encoding="utf-8") as handle:
            return json.load(handle)

    return {"sources": [], "filters": [], "sensors": []}


def load_component_database(filepath="component_database.json"):
    path = Path(filepath)
    if not path.exists():
        return default_component_database()

    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    defaults = default_component_database()
    for key in ("sources", "filters", "sensors"):
        data.setdefault(key, defaults[key])
    return data


def evaluate_source(source, wavelengths_nm):
    source_type = source.get("type", "flat")
    if source_type == "blackbody":
        return bounded_blackbody_spectrum(
            wavelengths_nm,
            float(source.get("temp_k", 3000)),
            min_nm=source.get("min_nm"),
            max_nm=source.get("max_nm"),
        )
    if source_type == "gaussian":
        return gaussian_bandpass(
            wavelengths_nm,
            float(source["center_nm"]),
            float(source["fwhm_nm"]),
            float(source.get("peak", 1.0)),
        )
    return np.full_like(wavelengths_nm, float(source.get("level", 1.0)), dtype=float)


def evaluate_filter(filter_def, wavelengths_nm):
    return gaussian_bandpass(
        wavelengths_nm,
        float(filter_def["center_nm"]),
        float(filter_def["fwhm_nm"]),
        float(filter_def.get("peak", 1.0)),
    )


def evaluate_sensor(sensor, wavelengths_nm):
    sensor_type = sensor.get("type", "flat")
    if sensor_type == "piecewise":
        return piecewise_response(wavelengths_nm, sensor["points"])
    return np.full_like(wavelengths_nm, float(sensor.get("level", 1.0)), dtype=float)
