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
    return {
        "sources": [
            {"name": "Blackbody Halogen 3000K", "type": "blackbody", "temp_k": 3000},
            {"name": "Blackbody Halogen 2800K", "type": "blackbody", "temp_k": 2800},
            {"name": "MTE6114W-WRC LED 1460nm", "type": "gaussian", "center_nm": 1460, "fwhm_nm": 103, "peak": 1.0},
            {"name": "HPIR104 Thermal Emitter", "type": "blackbody", "temp_k": 903.15, "min_nm": 2000, "max_nm": 11000},
            {"name": "Flat Ideal Source", "type": "flat", "level": 1.0},
        ],
        "filters": [
            {"name": "Thorlabs FL1064-10", "center_nm": 1064, "fwhm_nm": 10, "peak": 1.0},
            {"name": "Thorlabs FB1200-10", "center_nm": 1200, "fwhm_nm": 10, "peak": 1.0},
            {"name": "Thorlabs FB1400-12", "center_nm": 1400, "fwhm_nm": 12, "peak": 1.0},
            {"name": "Thorlabs FBH1450-12", "center_nm": 1450, "fwhm_nm": 12, "peak": 1.0},
            {"name": "Syron BP1535-6", "center_nm": 1535, "fwhm_nm": 6, "peak": 1.0},
            {"name": "Thorlabs FB1625-12", "center_nm": 1625, "fwhm_nm": 12, "peak": 1.0},
            {"name": "PE 1730 nm", "center_nm": 1730, "fwhm_nm": 18, "peak": 1.0},
            {"name": "PE 1750 nm", "center_nm": 1750, "fwhm_nm": 12, "peak": 1.0},
            {"name": "Water 1940 nm", "center_nm": 1940, "fwhm_nm": 30, "peak": 1.0},
            {"name": "EVOH 2012 nm", "center_nm": 2012, "fwhm_nm": 18, "peak": 1.0},
            {"name": "Thorlabs FB2100-12", "center_nm": 2100, "fwhm_nm": 12, "peak": 1.0},
            {"name": "Reference 2200 nm", "center_nm": 2200, "fwhm_nm": 20, "peak": 1.0},
            {"name": "PE 2310 nm", "center_nm": 2310, "fwhm_nm": 20, "peak": 1.0},
            {"name": "Andover 3000.0 / 100.0 - 69350-B", "center_nm": 3000, "fwhm_nm": 100, "peak": 1.0},
            {"name": "PIRF-INBP3400 Type 1", "center_nm": 3400, "fwhm_nm": 140, "peak": 1.0},
        ],
        "sensors": [
            {
                "name": "InGaAs",
                "type": "piecewise",
                "points": [[1500, 0.6], [2000, 0.9], [2300, 1.0], [2500, 0.8], [2600, 0.0]],
            },
            {
                "name": "InAsSb (2.7-5.3um)",
                "type": "piecewise",
                "points": [[2600, 0.0], [2700, 0.75], [3300, 1.0], [5000, 0.95], [5300, 0.5], [5400, 0.0]],
            },
            {
                "name": "MCT (MIR)",
                "type": "piecewise",
                "points": [[2500, 0.0], [2800, 0.8], [3200, 1.0], [3800, 0.9], [4100, 0.0]],
            },
            {"name": "Ideal (Flat)", "type": "flat", "level": 1.0},
        ],
    }


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
