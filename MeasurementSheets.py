"""Read measurement tables (.xlsx/.csv) and turn selected cell ranges into data points.

No UI dependencies: the viewer in ``BroadbandSheetUI.py`` only displays what
``read_tabular_file`` returns and calls ``measurement_points`` with the cell
range the user picked.  Cells are kept as display strings and are only
converted to numbers here, where every rejection can be explained.
"""

from __future__ import annotations

import csv
import os

# A sheet is shown, not analysed, so a generous cap keeps the viewer responsive
# without hiding realistic measurement tables.
MAX_DISPLAY_ROWS = 300
MAX_DISPLAY_COLUMNS = 30


def column_label(index):
    """0 -> 'A', 25 -> 'Z', 26 -> 'AA'."""
    if index < 0:
        raise ValueError("Column index must be zero or more.")
    label = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        label = chr(ord("A") + remainder) + label
    return label


def cell_label(row, column):
    """(0, 1) -> 'B1'."""
    return f"{column_label(column)}{row + 1}"


def range_label(top, left, bottom, right):
    """(1, 1, 4, 1) -> 'B2:B5'; a single cell is shown on its own."""
    first, last = cell_label(top, left), cell_label(bottom, right)
    return first if first == last else f"{first}:{last}"


def cells_in_range(grid, top, left, bottom, right):
    """Flat row-major list of display strings, padding short rows with ''."""
    values = []
    for row_index in range(top, bottom + 1):
        row = grid[row_index] if row_index < len(grid) else []
        for column_index in range(left, right + 1):
            values.append(row[column_index] if column_index < len(row) else "")
    return values


def _cell_text(value):
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        rows = [[cell.strip() for cell in row] for row in csv.reader(handle, dialect)]
    return {os.path.basename(path): rows}


def _read_workbook(path):
    try:
        from openpyxl import load_workbook
    except ImportError as exc:
        raise RuntimeError("Reading .xlsx files requires openpyxl. Install the project requirements first.") from exc
    # data_only: formula cells give their last saved result, not the formula text.
    workbook = load_workbook(path, data_only=True, read_only=True)
    sheets = {}
    try:
        for worksheet in workbook.worksheets:
            rows = []
            for row in worksheet.iter_rows(max_row=MAX_DISPLAY_ROWS, max_col=MAX_DISPLAY_COLUMNS,
                                           values_only=True):
                rows.append([_cell_text(value) for value in row])
            sheets[worksheet.title] = rows
    finally:
        workbook.close()
    return sheets


def read_tabular_file(path):
    """Read a spreadsheet into {sheet name: rows of display strings}.

    Returns a mapping with ``sheets``, the display cap that was applied and
    whether anything was cut off, so the viewer can say so.
    """
    extension = os.path.splitext(path)[1].lower()
    if extension in (".csv", ".txt", ".tsv"):
        sheets = _read_csv(path)
    elif extension in (".xlsx", ".xlsm", ".xltx"):
        sheets = _read_workbook(path)
    else:
        raise ValueError(f"{os.path.basename(path)}: open a .xlsx or .csv file "
                         f"(.xls must be re-saved as .xlsx first).")
    if not sheets or all(not rows for rows in sheets.values()):
        raise ValueError(f"{os.path.basename(path)} has no readable rows.")

    truncated = False
    for name, rows in sheets.items():
        if len(rows) > MAX_DISPLAY_ROWS:
            rows[:] = rows[:MAX_DISPLAY_ROWS]
            truncated = True
        for row in rows:
            if len(row) > MAX_DISPLAY_COLUMNS:
                del row[MAX_DISPLAY_COLUMNS:]
                truncated = True
        # Trailing blank rows and cells only pad the viewer; openpyxl pads every
        # row out to max_col, which would show a 2-column sheet as 30 columns.
        for row in rows:
            while row and not row[-1]:
                row.pop()
        while rows and not any(cell for cell in rows[-1]):
            rows.pop()
    return {"path": path, "name": os.path.basename(path), "sheets": sheets,
            "truncated": truncated, "max_rows": MAX_DISPLAY_ROWS, "max_columns": MAX_DISPLAY_COLUMNS}


def _number(text):
    """Parse one cell, tolerating thousands separators and a trailing unit."""
    cleaned = str(text).strip().replace(",", "")
    if not cleaned:
        return None
    for unit in ("mv", "mV", "µm", "um", "layers", "layer"):
        if cleaned.lower().endswith(unit.lower()):
            cleaned = cleaned[: -len(unit)].strip()
            break
    try:
        return float(cleaned)
    except ValueError:
        return None


def measurement_points(thickness_cells, voltage_cells, thickness_label="thickness", voltage_label="voltage"):
    """Pair two selected ranges into measurement points, explaining every rejection.

    * Pairs where BOTH cells are non-numeric are dropped as headers or blanks.
    * A pair where only ONE cell is numeric is an error: the two selections are
      misaligned, and guessing which row to drop would silently shift the data.
    * A zero-thickness row is returned separately as the no-film reading.
    """
    if len(thickness_cells) != len(voltage_cells):
        raise ValueError(f"The two selections must cover the same number of cells "
                         f"({len(thickness_cells)} {thickness_label} vs {len(voltage_cells)} {voltage_label}).")
    if not thickness_cells:
        raise ValueError("Select the cells holding the measurements first.")

    points, skipped, no_film = [], 0, None
    for index, (raw_thickness, raw_voltage) in enumerate(zip(thickness_cells, voltage_cells)):
        thickness, voltage = _number(raw_thickness), _number(raw_voltage)
        if thickness is None and voltage is None:
            skipped += 1
            continue
        if thickness is None or voltage is None:
            missing = thickness_label if thickness is None else voltage_label
            text = raw_thickness if thickness is None else raw_voltage
            raise ValueError(f"Row {index + 1} of the selection has a {missing} of {text!r}, which is not a number, "
                             f"while the other column has one. Re-select so the two ranges line up.")
        if thickness < 0:
            raise ValueError(f"Row {index + 1} of the selection has a negative {thickness_label} ({thickness:g}).")
        if voltage < 0:
            raise ValueError(f"Row {index + 1} of the selection has a negative {voltage_label} ({voltage:g} mV).")
        if thickness == 0:
            if no_film is not None and no_film != voltage:
                raise ValueError(f"The selection has two different zero-thickness readings "
                                 f"({no_film:g} and {voltage:g} mV).")
            no_film = voltage
            continue
        points.append((thickness, voltage))

    duplicates = {value for value, _ in points if [v for v, _ in points].count(value) > 1}
    if duplicates:
        raise ValueError(f"The selection repeats the same {thickness_label} "
                         f"({', '.join(f'{value:g}' for value in sorted(duplicates))}). "
                         f"Average the repeats before importing them.")
    if not points:
        raise ValueError("No numeric measurement pairs were found in the selection.")
    points.sort()
    return {"points": points, "skipped_rows": skipped, "no_film_voltage_mv": no_film}
