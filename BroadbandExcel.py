"""Human-readable Excel export for broadband detector analysis."""

from __future__ import annotations

from pathlib import Path

# Per-layer spectral columns are written for at most this many layer counts so a
# 1000-point sweep does not produce an unusably wide sheet.
MAX_SPECTRAL_LAYER_COLUMNS = 20


def _number(value):
    """openpyxl cannot write inf/nan; show them as blank cells instead."""
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return value
    return value if value == value and abs(value) != float("inf") else None


def export_broadband_workbook(path, analysis, metadata, experimental_rows=None):
    try:
        from openpyxl import Workbook
        from openpyxl.chart import BarChart, Reference, ScatterChart, Series
        from openpyxl.styles import Alignment, Font, PatternFill
        from openpyxl.utils import get_column_letter
    except ImportError as exc:
        raise RuntimeError("Excel export requires openpyxl. Install project requirements first.") from exc

    path = Path(path)
    if path.suffix.lower() != ".xlsx":
        path = path.with_suffix(".xlsx")
    wb = Workbook()
    header_fill = PatternFill("solid", fgColor="1F4E78")
    section_fill = PatternFill("solid", fgColor="D9E2F3")
    header_font = Font(color="FFFFFF", bold=True)
    bold = Font(bold=True)
    rows = analysis["rows"]
    n_rows = len(rows)

    def style_header(sheet, row_index):
        for cell in sheet[row_index]:
            if cell.value is not None:
                cell.fill, cell.font = header_fill, header_font
                cell.alignment = Alignment(wrap_text=True, vertical="center")

    def section(sheet, title):
        sheet.append([title])
        cell = sheet.cell(row=sheet.max_row, column=1)
        cell.font, cell.fill = Font(bold=True, size=12), section_fill

    # ---------------------------------------------------------------- Summary
    ws = wb.active
    ws.title = "Summary"
    ws.append(["Broadband Spectral Detector Analysis"])
    ws["A1"].font = Font(size=16, bold=True)
    ws.append(["Absorbance spectrum → Beer-Lambert per wavelength → detector-weighted integral → predicted voltage"])
    ws.append([])

    section(ws, "Experiment settings")
    for label, key in [
        ("Material", "material"), ("Spectrum file", "spectrum_source"), ("Detector", "detector"),
        ("Detector weighting", "weighting_mode"), ("Integration range (nm)", "wavelength_range"),
        ("No-film voltage V0 (mV)", "no_film_voltage_mv"), ("Simulation input", "simulation_input"),
        ("Reference thickness", "reference_thickness"), ("Absorbance convention", "absorbance_mode"),
        ("Baseline correction", "baseline_correction"), ("Baseline offset removed (A)", "baseline_offset"),
        ("Negative absorbance clamped", "clamp_negative"), ("Fresnel reflection loss", "fresnel"),
        ("Interface transmission per layer", "layer_interface_transmission"),
        ("Assumptions / notes", "notes"),
    ]:
        value = metadata.get(key, "")
        if isinstance(value, bool):
            value = "Yes" if value else "No"
        ws.append([label, value if isinstance(value, str) else _number(value)])
        ws.cell(row=ws.max_row, column=1).font = bold
        ws.cell(row=ws.max_row, column=2).alignment = Alignment(wrap_text=True, vertical="top")
    ws.append([])

    if metadata.get("inspected_scale") is not None:
        section(ws, f"Headline result ({metadata.get('simulation_input', 'Layer count')} = {metadata['inspected_scale']:g})")
        for label, key in [("Measured (mV)", "inspected_measured_mv"), ("Predicted from spectral model (mV)", "inspected_predicted_mv"),
                           ("Error (mV)", "inspected_error_mv"), ("Error (%)", "inspected_error_percent")]:
            ws.append([label, _number(metadata.get(key))])
            ws.cell(row=ws.max_row, column=1).font = bold
        ws.append([])

    section(ws, "Measured vs predicted")
    ws.append(["Layer Count", "Measured Voltage (mV)", "Predicted Voltage (mV)", "Error (mV)",
               "Absolute Error (mV)", "Percent Error (%)", "Absolute Percent Error (%)",
               "Measured Transmission V/V0", "Predicted Transmission V/V0", "Spectral Transmission",
               "Interface Transmission", "ln(V0/V) Measured", "ln(V0/V) Predicted"])
    header_row = ws.max_row
    style_header(ws, header_row)
    for row in rows:
        ws.append([_number(row["thickness_scale"]), _number(row["measured_voltage_mv"]),
                   _number(row["predicted_voltage_mv"]), _number(row["error_mv"]),
                   _number(row["absolute_error_mv"]), _number(row["percent_error"]),
                   _number(row["absolute_percent_error"]), _number(row["measured_transmission"]),
                   _number(row["predicted_transmission"]), _number(row.get("spectral_transmission")),
                   _number(row.get("interface_transmission")), _number(row["measured_attenuation_natural"]),
                   _number(row["predicted_attenuation_natural"])])
    first_data, last_data = header_row + 1, header_row + n_rows
    for r in range(first_data, last_data + 1):
        for c in range(2, 8):
            ws.cell(row=r, column=c).number_format = "0.0"
        for c in range(8, 14):
            ws.cell(row=r, column=c).number_format = "0.0000"
    ws.append([])
    section(ws, "Model agreement metrics")
    metrics = analysis["metrics"]
    for label, key in [("MAE (mV)", "mae_mv"), ("RMSE (mV)", "rmse_mv"),
                       ("MAPE (%)", "mape_percent"), ("R² (measured vs predicted voltage)", "r_squared")]:
        ws.append([label, _number(metrics.get(key))])
        ws.cell(row=ws.max_row, column=1).font = bold

    # ------------------------------------------------- Spectral calculation
    spectral = analysis["spectral"]
    spectra_by_scale = analysis.get("spectra_by_scale", {spectral["thickness_scale"]: spectral})
    layer_scales = list(spectra_by_scale)[:MAX_SPECTRAL_LAYER_COLUMNS]
    wavelength = spectral["wavelength_nm"]
    n = len(wavelength)
    raw = spectral.get("raw_absorbance", spectral["absorbance"])
    source = spectral.get("source")
    response = spectral.get("responsivity")

    spec_ws = wb.create_sheet("Spectral Calculation")
    base_headers = ["Wavelength (nm)", "Raw Absorbance", "Absorbance Used (1 layer)", "Source S(λ)",
                    "Detector Responsivity R(λ)", "Weight W(λ)"]
    layer_headers = []
    for scale in layer_scales:
        layer_headers += [f"T(λ) @ {scale:g}", f"W·T @ {scale:g}"]
    spec_ws.append(base_headers + layer_headers)
    style_header(spec_ws, 1)
    for i in range(n):
        values = [float(wavelength[i]), float(raw[i]), float(spectral["absorbance"][i]),
                  None if source is None else float(source[i]),
                  None if response is None else float(response[i]),
                  float(spectral["weight"][i])]
        for scale in layer_scales:
            result = spectra_by_scale[scale]
            interface = result.get("interface_transmission", 1.0)
            values += [_number(result["transmission"][i] * interface),
                       _number(result["weighted_transmitted"][i] * interface)]
        spec_ws.append(values)

    # ------------------------------------------------------ Layer predictions
    layer_ws = wb.create_sheet("Layer Predictions")
    layer_ws.append(["Layer Count", "∫W(λ)dλ", "∫W(λ)T(λ)dλ", "Spectral Transmission",
                     "Interface Transmission", "Effective Transmission", "Effective Attenuation ln(V0/V)",
                     "Effective Absorbance -log10(V/V0)", "Predicted Voltage (mV)"])
    style_header(layer_ws, 1)
    try:
        import numpy as np
        trapezoid = getattr(np, "trapezoid", None) or np.trapz
    except ImportError:  # pragma: no cover - numpy is a hard dependency elsewhere.
        trapezoid = None
    for row in rows:
        result = spectra_by_scale.get(row["thickness_scale"])
        integral_w = integral_wt = None
        if result is not None and trapezoid is not None:
            integral_w = float(trapezoid(result["weight"], result["wavelength_nm"]))
            integral_wt = float(trapezoid(result["weighted_transmitted"], result["wavelength_nm"]))
        layer_ws.append([_number(row["thickness_scale"]), integral_w, integral_wt,
                         _number(row.get("spectral_transmission")), _number(row.get("interface_transmission")),
                         _number(row["predicted_transmission"]), _number(row["predicted_attenuation_natural"]),
                         _number(row["predicted_absorbance_base10"]), _number(row["predicted_voltage_mv"])])

    # ------------------------------------------------------ Experimental data
    exp_ws = wb.create_sheet("Experimental Data")
    exp_ws.append(["Material", "Detector", "Layer Count", "Measured Voltage (mV)", "No-film Voltage (mV)",
                   "Measured V/V0"])
    style_header(exp_ws, 1)
    for row in experimental_rows or []:
        v0 = row.get("no_film_voltage_mv")
        measured = row.get("measured_voltage_mv")
        exp_ws.append([row.get("material"), row.get("detector"), row.get("layer_count"), measured, v0,
                       measured / v0 if measured is not None and v0 else None])

    # ----------------------------------------------------------------- Charts
    chart_ws = wb.create_sheet("Charts")

    def scatter(title, x_title, y_title):
        chart = ScatterChart()
        chart.title, chart.style = title, 13
        chart.x_axis.title, chart.y_axis.title = x_title, y_title
        chart.height, chart.width = 9, 16
        # Excel hides axes on new scatter charts unless delete is explicitly off.
        chart.x_axis.delete = False
        chart.y_axis.delete = False
        return chart

    def add_series(chart, sheet, x_col, y_col, first, last, title, markers=False, line=True):
        series = Series(Reference(sheet, min_col=y_col, min_row=first, max_row=last),
                        Reference(sheet, min_col=x_col, min_row=first, max_row=last), title=title)
        series.smooth = False
        if markers:
            series.marker.symbol, series.marker.size = "circle", 7
        else:
            series.marker.symbol = "none"
        if not line:
            series.graphicalProperties.line.noFill = True
        chart.series.append(series)
        return series

    compare_chart = scatter("Measured vs Predicted Voltage", "Layer count", "Detector voltage (mV)")
    add_series(compare_chart, ws, 1, 3, first_data, last_data, "Predicted (spectral model)", markers=True)
    add_series(compare_chart, ws, 1, 2, first_data, last_data, "Measured", markers=True, line=False)
    chart_ws.add_chart(compare_chart, "A1")

    parity_chart = scatter("Parity: Predicted vs Measured", "Measured voltage (mV)", "Predicted voltage (mV)")
    add_series(parity_chart, ws, 2, 3, first_data, last_data, "Layer counts", markers=True, line=False)
    chart_ws.add_chart(parity_chart, "K1")

    error_chart = BarChart()
    error_chart.title, error_chart.height, error_chart.width = "Prediction Error (%)", 9, 16
    error_chart.y_axis.title, error_chart.x_axis.title = "Error (%)  + = over-predicts", "Layer count"
    error_chart.x_axis.delete = False
    error_chart.y_axis.delete = False
    error_chart.add_data(Reference(ws, min_col=6, min_row=header_row, max_row=last_data), titles_from_data=True)
    error_chart.set_categories(Reference(ws, min_col=1, min_row=first_data, max_row=last_data))
    chart_ws.add_chart(error_chart, "A20")

    attenuation_chart = scatter("Effective Attenuation vs Layers", "Layer count", "ln(V0/V)")
    add_series(attenuation_chart, ws, 1, 13, first_data, last_data, "Predicted", markers=True)
    add_series(attenuation_chart, ws, 1, 12, first_data, last_data, "Measured", markers=True, line=False)
    chart_ws.add_chart(attenuation_chart, "K20")

    absorbance_chart = scatter("Absorbance Spectrum", "Wavelength (nm)", "Absorbance")
    add_series(absorbance_chart, spec_ws, 1, 2, 2, n + 1, "Raw")
    add_series(absorbance_chart, spec_ws, 1, 3, 2, n + 1, "Used by model")
    chart_ws.add_chart(absorbance_chart, "A39")

    transmission_chart = scatter("Transmission Spectrum per Layer Count", "Wavelength (nm)", "Transmission")
    for index, scale in enumerate(layer_scales):
        add_series(transmission_chart, spec_ws, 1, 7 + 2 * index, 2, n + 1, f"{scale:g}")
    chart_ws.add_chart(transmission_chart, "K39")

    selected = analysis.get("selected_scale", layer_scales[0])
    selected_index = layer_scales.index(selected) if selected in layer_scales else 0
    contribution_chart = scatter(f"Wavelengths Contributing to Detector ({layer_scales[selected_index]:g})",
                                 "Wavelength (nm)", "Weighted contribution")
    add_series(contribution_chart, spec_ws, 1, 6, 2, n + 1, "No film W(λ)")
    add_series(contribution_chart, spec_ws, 1, 8 + 2 * selected_index, 2, n + 1, "With film W·T")
    chart_ws.add_chart(contribution_chart, "A58")

    # ----------------------------------------------------------------- Layout
    ws.column_dimensions["A"].width = 36
    ws.column_dimensions["B"].width = 44
    for column_index in range(3, 14):
        ws.column_dimensions[get_column_letter(column_index)].width = 16
    for sheet in (spec_ws, layer_ws, exp_ws):
        sheet.freeze_panes = "B2" if sheet is spec_ws else "A2"
        for column_index in range(1, sheet.max_column + 1):
            header = sheet.cell(row=1, column=column_index).value
            sheet.column_dimensions[get_column_letter(column_index)].width = max(14, min(32, len(str(header or "")) + 4))
    wb.save(path)
    return str(path)
