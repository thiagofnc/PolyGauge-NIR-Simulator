"""File loaders for broadband analysis: material spectra, weighting curves and presets.

No UI dependencies.  Every loader returns wavelengths in nanometres and states
what it assumed, and none of them silently rescales or clips the data.
"""

from __future__ import annotations

import csv
import json
import os

import numpy as np

from MeasuredData import discover_measured_samples, load_log_spectrum, wavenumber_to_nm

MATERIAL_METADATA_FILE = "material_spectra.json"
DETECTOR_DATA_FILE = "experimental_detector_data.json"
COMPONENT_DATABASE_FILE = "component_database.json"


def _sorted_xy(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0)
    order = np.argsort(x[valid], kind="stable")
    return x[valid][order], y[valid][order]


def read_two_column_csv(path):
    """Read the first two numeric columns.  Returns (x, y, header_or_None)."""
    xs, ys, header = [], [], None
    with open(path, newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t ")
        except csv.Error:
            dialect = csv.excel
        for row in csv.reader(handle, dialect):
            cells = [cell.strip() for cell in row if cell.strip()]
            if len(cells) < 2:
                continue
            try:
                xs.append(float(cells[0]))
                ys.append(float(cells[1]))
            except ValueError:
                if not xs and header is None:
                    header = cells[:2]
                    continue
                raise ValueError(f"{os.path.basename(path)}: non-numeric row {cells[:2]}")
    if len(xs) < 2:
        raise ValueError(f"{os.path.basename(path)} needs at least two numeric rows.")
    return np.asarray(xs), np.asarray(ys), header


def wavelength_unit_from_header(header):
    """'nm', 'um' or 'cm-1' from a column header; 'nm' when there is no header."""
    if not header:
        return "nm", "No header row; wavelength assumed to be nm."
    text = header[0].lower().replace("μ", "µ")
    if "cm-1" in text or "cm^-1" in text or "1/cm" in text or "wavenumber" in text:
        return "cm-1", f"Header '{header[0]}' read as wavenumber (cm⁻¹)."
    if "µm" in text or "um)" in text or text.endswith("um") or "micron" in text:
        return "um", f"Header '{header[0]}' read as µm."
    if "nm" in text:
        return "nm", f"Header '{header[0]}' read as nm."
    return "nm", f"Header '{header[0]}' has no unit; wavelength assumed to be nm."


def to_nm(x, unit):
    x = np.asarray(x, dtype=float)
    if unit == "nm":
        return x
    if unit == "um":
        return x * 1000.0
    if unit == "cm-1":
        return wavenumber_to_nm(x)
    raise ValueError(f"Unknown wavelength unit {unit!r}.")


def load_spectral_curve_csv(path, kind):
    """Load a detector responsivity, source spectrum or optical-term CSV.

    ``kind`` is 'responsivity', 'source' or 'transmission'.  Columns:
    wavelength, value.  Negative values are rejected instead of clipped, and a
    transmission term must be a 0-1 fraction.
    """
    x, y, header = read_two_column_csv(path)
    unit, unit_note = wavelength_unit_from_header(header)
    wavelength, values = _sorted_xy(to_nm(x, unit), y)
    if np.any(values < 0):
        raise ValueError(f"{os.path.basename(path)} contains negative {kind} values; clean the file first.")
    if kind == "transmission" and np.max(values) > 1.0 + 1e-6:
        raise ValueError(f"{os.path.basename(path)}: transmission terms must be fractions 0-1, not percent.")
    if np.max(values) <= 0:
        raise ValueError(f"{os.path.basename(path)} has no positive {kind} values.")
    return {
        "wavelength_nm": wavelength, "values": values, "path": path, "kind": kind,
        "note": f"{os.path.basename(path)}: {wavelength.size} points, {wavelength[0]:.0f}-{wavelength[-1]:.0f} nm. "
                f"{unit_note} Zero outside this range.",
    }


def _transmittance_to_absorbance(values, quantity, name):
    values = np.asarray(values, dtype=float)
    if quantity == "transmittance_fraction":
        if np.nanmax(values) > 1.5:
            raise ValueError(f"{name}: declared as 0-1 transmittance but values exceed 1.5 (percent?).")
        fraction = values
    elif quantity == "transmittance_percent":
        fraction = values / 100.0
    else:
        raise ValueError(f"{name}: unknown quantity {quantity!r}.")
    if np.any(fraction <= 0):
        raise ValueError(f"{name}: transmittance must be positive to convert to absorbance.")
    return -np.log10(fraction)


def load_material_spectrum(base_dir, name, entry):
    """Load one material described in material_spectra.json."""
    path = os.path.join(base_dir, entry["file"])
    quantity = entry["quantity"]
    if entry["format"] == "jcamp_log":
        wavenumber, values, header = load_log_spectrum(path)
        yunits = header.get("YUNITS", "").strip().lower()
        if quantity.startswith("absorbance") and yunits and not yunits.startswith("a"):
            raise ValueError(f"{name}: declared absorbance but file YUNITS is {header.get('YUNITS')!r}.")
        wavelength, values = _sorted_xy(wavenumber_to_nm(wavenumber), values)
    elif entry["format"] == "csv_wavelength_nm":
        x, values, _header = read_two_column_csv(path)
        wavelength, values = _sorted_xy(x, values)
    else:
        raise ValueError(f"{name}: unknown format {entry['format']!r}.")

    if quantity == "absorbance_base10":
        convention, absorbance = "base10", values
    elif quantity == "absorbance_natural":
        convention, absorbance = "natural", values
    else:
        convention, absorbance = "base10", _transmittance_to_absorbance(values, quantity, name)
    reference = entry.get("reference_thickness_um")
    return {
        "name": name, "wavelength_nm": wavelength, "absorbance": absorbance, "convention": convention,
        "reference_thickness_um": None if reference is None else float(reference),
        "source": entry["file"], "note": entry.get("note", ""),
    }


def load_material_library(base_dir, include_discovered=True):
    """Named materials from material_spectra.json plus every discovered FTIR absorbance log."""
    materials = {}
    metadata_path = os.path.join(base_dir, MATERIAL_METADATA_FILE)
    if os.path.exists(metadata_path):
        with open(metadata_path, encoding="utf-8") as handle:
            for name, entry in json.load(handle).items():
                if os.path.exists(os.path.join(base_dir, entry["file"])):
                    materials[name] = load_material_spectrum(base_dir, name, entry)
    if not include_discovered:
        return materials

    known = {os.path.normcase(os.path.abspath(os.path.join(base_dir, m["source"]))) for m in materials.values()}
    samples, _headers = discover_measured_samples(os.path.join(base_dir, "logs_full_range"),
                                                  os.path.join(base_dir, "sample_references.xlsx"))
    for sample in samples:
        path = sample["paths"].get("absorbance")
        if not path or os.path.normcase(os.path.abspath(path)) in known:
            continue
        wavenumber, values, header = load_log_spectrum(path)
        if not header.get("YUNITS", "abs").strip().lower().startswith("a"):
            continue
        wavelength, values = _sorted_xy(wavenumber_to_nm(wavenumber), values)
        name = f"Measured: {sample['label']}"
        materials[name] = {
            "name": name, "wavelength_nm": wavelength, "absorbance": values, "convention": "base10",
            "reference_thickness_um": None, "source": os.path.relpath(path, base_dir),
            "note": "Discovered FTIR absorbance log (YUNITS=Abs, base-10). Reference thickness not recorded.",
        }
    return materials


def load_detector_definitions(base_dir):
    with open(os.path.join(base_dir, DETECTOR_DATA_FILE), encoding="utf-8") as handle:
        return json.load(handle)


def load_source_presets(base_dir):
    """Blackbody emitters from component_database.json as {label: preset}."""
    presets = {}
    path = os.path.join(base_dir, COMPONENT_DATABASE_FILE)
    if not os.path.exists(path):
        return presets
    with open(path, encoding="utf-8") as handle:
        database = json.load(handle)
    for source in database.get("sources", []):
        if source.get("type") != "blackbody" or not source.get("temp_k"):
            continue
        label = f"{source['name']} ({source['temp_k']:g} K"
        if source.get("min_nm") and source.get("max_nm"):
            label += f", {source['min_nm']:g}-{source['max_nm']:g} nm"
        presets[label + ")"] = {
            "name": source["name"], "temperature_k": float(source["temp_k"]),
            "min_nm": source.get("min_nm"), "max_nm": source.get("max_nm"),
        }
    return presets
