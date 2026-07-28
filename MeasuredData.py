"""Loaders for the measured spectrometer logs and the sample reference sheet.

The logs in ``logs_full_range/`` are JCAMP-style text exports from the FTIR: four
``##KEY=VALUE`` header lines followed by ``<wavenumber cm^-1>\t<value>`` rows.
Files ending in ``_abs`` carry absorbance, ``_trans`` carry percent
transmittance.  ``sample_references.xlsx`` describes what each physical sample
is (layers, thickness, filler, composition).

The xlsx is parsed with the standard library (zipfile + ElementTree) so the app
does not gain an openpyxl dependency just to read one 16-row table.
"""

import os
import re
import zipfile
from xml.etree import ElementTree as ET

import numpy as np

LOG_DIR = "logs_full_range"
REFERENCE_XLSX = "sample_references.xlsx"

ABS_TOKENS = ("abs", "absorbance")
TRANS_TOKENS = ("trans", "transmittance", "transmission")

_SHEET_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


# --- Reference sheet -------------------------------------------------------

def _cell_column(ref):
    """'C12' -> 'C'."""
    return "".join(ch for ch in ref if ch.isalpha())


def read_reference_table(path=REFERENCE_XLSX):
    """Returns (headers, rows_by_sample_number).

    ``headers`` is an ordered list of (column_letter, header_text) for the
    table's header row; ``rows_by_sample_number`` maps 1 -> {header: value}.
    Returns ([], {}) if the workbook is missing or unreadable.
    """
    if not os.path.exists(path):
        return [], {}

    try:
        with zipfile.ZipFile(path) as zf:
            shared = []
            if "xl/sharedStrings.xml" in zf.namelist():
                root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
                for si in root:
                    shared.append("".join(t.text or "" for t in si.iter(_SHEET_NS + "t")))
            sheet = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
    except (zipfile.BadZipFile, KeyError, ET.ParseError):
        return [], {}

    grid = {}
    for row in sheet.iter(_SHEET_NS + "row"):
        cells = {}
        for cell in row.iter(_SHEET_NS + "c"):
            value_node = cell.find(_SHEET_NS + "v")
            if value_node is None or value_node.text is None:
                continue
            if cell.get("t") == "s":
                text = shared[int(value_node.text)]
            elif cell.get("t") == "inlineStr":
                text = "".join(t.text or "" for t in cell.iter(_SHEET_NS + "t"))
            else:
                text = value_node.text
            text = text.strip()
            if text:
                cells[_cell_column(cell.get("r", ""))] = text
        if cells:
            grid[int(row.get("r"))] = cells

    # The header row is the one whose first populated cell reads "Sample name".
    header_row = None
    for row_index in sorted(grid):
        if any(v.lower() == "sample name" for v in grid[row_index].values()):
            header_row = row_index
            break
    if header_row is None:
        return [], {}

    headers = [(col, grid[header_row][col]) for col in sorted(grid[header_row])]

    rows = {}
    for row_index in sorted(grid):
        if row_index <= header_row:
            continue
        cells = grid[row_index]
        name_col = headers[0][0]
        match = re.match(r"sample\s*(\d+)", cells.get(name_col, ""), re.IGNORECASE)
        if not match:
            continue
        rows[int(match.group(1))] = {
            header: cells.get(col, "") for col, header in headers
        }
    return headers, rows


# --- Log files -------------------------------------------------------------

def parse_log_name(filename):
    """Splits a log filename into (sample_number, variant, mode).

    ``mode`` is 'absorbance' or 'transmittance'; ``variant`` is whatever else
    the name carries ('no_film', 'label_down_new', ...) and is '' for a plain
    ``sampleN_abs.txt``.  Returns None when the name has no abs/trans marker.
    """
    stem = os.path.splitext(os.path.basename(filename))[0]
    mode = None
    rest = []
    for part in stem.split("_"):
        low = part.lower()
        if low in ABS_TOKENS:
            mode = "absorbance"
        elif low in TRANS_TOKENS:
            mode = "transmittance"
        else:
            rest.append(part)
    if mode is None or not rest:
        return None

    match = re.fullmatch(r"sample\s*(\d+)", rest[0], re.IGNORECASE)
    if not match:
        return None
    return int(match.group(1)), "_".join(rest[1:]), mode


def discover_measured_samples(log_dir=LOG_DIR, reference_path=REFERENCE_XLSX):
    """Builds the list of measured samples available for plotting.

    Each entry is a dict with the sample number, the variant suffix, both file
    paths (either may be None), the matching reference-sheet row and a display
    label.  Sorted by sample number, then variant.
    """
    headers, reference_rows = read_reference_table(reference_path)
    if not os.path.isdir(log_dir):
        return [], headers

    grouped = {}
    for filename in sorted(os.listdir(log_dir)):
        if not filename.lower().endswith(".txt"):
            continue
        parsed = parse_log_name(filename)
        if parsed is None:
            continue
        number, variant, mode = parsed
        entry = grouped.setdefault((number, variant), {
            "number": number,
            "variant": variant,
            "paths": {"absorbance": None, "transmittance": None},
        })
        entry["paths"][mode] = os.path.join(log_dir, filename)

    samples = []
    for (number, variant) in sorted(grouped, key=lambda k: (k[0], k[1])):
        entry = grouped[(number, variant)]
        reference = reference_rows.get(number, {})
        entry["reference"] = reference
        entry["key"] = f"sample{number}" + (f"_{variant}" if variant else "")
        entry["label"] = _build_label(number, variant, reference, headers)
        samples.append(entry)
    return samples, headers


def _build_label(number, variant, reference, headers):
    label = f"Sample {number}"
    if variant:
        label += f" ({variant.replace('_', ' ')})"
    # The last reference column is the shorthand composition name.
    if headers:
        composition = reference.get(headers[-1][1], "")
        if composition:
            label += f" - {composition}"
    return label


def load_log_spectrum(path):
    """Reads a JCAMP-style log into (wavenumber_cm1, values, header_dict).

    Rows are returned sorted by ascending wavenumber; malformed lines are
    skipped.
    """
    header = {}
    x_values = []
    y_values = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("##"):
                key, _, value = line[2:].partition("=")
                header[key.strip().upper()] = value.strip()
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue
            try:
                x_values.append(float(parts[0]))
                y_values.append(float(parts[1]))
            except ValueError:
                continue

    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    if x.size:
        order = np.argsort(x)
        x, y = x[order], y[order]
    return x, y, header


def wavenumber_to_nm(wavenumber_cm1):
    """cm^-1 -> nm, guarding against a zero entry."""
    wavenumber_cm1 = np.asarray(wavenumber_cm1, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        nm = 1.0e7 / wavenumber_cm1
    return nm
