"""Tests for reading measurement tables and turning cell ranges into data points."""

import os
import tempfile
import unittest

from MeasurementSheets import (cell_label, cells_in_range, column_label, measurement_points, range_label,
                               read_tabular_file)

GRID = [["Layers", "Voltage (mV)", "note"],
        ["1", "456", "run A"],
        ["2", "418", ""],
        ["3", "384", ""],
        ["4", "340", ""]]


def write_workbook(directory, rows, sheet_names=("Measurements",)):
    from openpyxl import Workbook
    workbook = Workbook()
    workbook.remove(workbook.active)
    for name in sheet_names:
        sheet = workbook.create_sheet(name)
        for row in rows:
            sheet.append(row)
    path = os.path.join(directory, "measurements.xlsx")
    workbook.save(path)
    return path


class LabelTests(unittest.TestCase):
    def test_column_labels_go_past_z(self):
        self.assertEqual([column_label(i) for i in (0, 1, 25, 26, 27, 51, 52)],
                         ["A", "B", "Z", "AA", "AB", "AZ", "BA"])
        with self.assertRaises(ValueError):
            column_label(-1)

    def test_cell_and_range_labels(self):
        self.assertEqual(cell_label(0, 0), "A1")
        self.assertEqual(cell_label(4, 1), "B5")
        self.assertEqual(range_label(1, 1, 4, 1), "B2:B5")
        self.assertEqual(range_label(2, 0, 2, 0), "A3")

    def test_cells_in_range_is_row_major_and_pads_short_rows(self):
        grid = [["a", "b"], ["c"]]
        self.assertEqual(cells_in_range(grid, 0, 0, 1, 1), ["a", "b", "c", ""])


class ReadTests(unittest.TestCase):
    def test_reads_every_sheet_of_a_workbook(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_workbook(tmp, GRID, sheet_names=("First", "Second"))
            book = read_tabular_file(path)
            self.assertEqual(list(book["sheets"]), ["First", "Second"])
            self.assertEqual(book["sheets"]["First"][1][:2], ["1", "456"])
            self.assertFalse(book["truncated"])
            self.assertEqual(book["name"], "measurements.xlsx")

    def test_numbers_keep_their_plain_form(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_workbook(tmp, [["Layers", "mV"], [1, 456.0], [2, 417.5]])
            rows = read_tabular_file(path)["sheets"]["Measurements"]
            self.assertEqual(rows[1], ["1", "456"])      # 456.0 is not shown as "456.0"
            self.assertEqual(rows[2], ["2", "417.5"])

    def test_reads_csv_as_a_single_sheet(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "m.csv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("Layers,Voltage\n1,456\n2,418\n")
            book = read_tabular_file(path)
            self.assertEqual(list(book["sheets"]), ["m.csv"])
            self.assertEqual(book["sheets"]["m.csv"][2], ["2", "418"])

    def test_trailing_blank_rows_are_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = write_workbook(tmp, [["Layers", "mV"], [1, 456], [None, None], [None, None]])
            self.assertEqual(len(read_tabular_file(path)["sheets"]["Measurements"]), 2)

    def test_unsupported_and_empty_files_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            legacy = os.path.join(tmp, "old.xls")
            open(legacy, "w").close()
            with self.assertRaises(ValueError) as caught:
                read_tabular_file(legacy)
            self.assertIn("re-saved as .xlsx", str(caught.exception))
            empty = os.path.join(tmp, "empty.csv")
            open(empty, "w").close()
            with self.assertRaises(ValueError):
                read_tabular_file(empty)


class MeasurementPointTests(unittest.TestCase):
    def test_pairs_two_selections_and_skips_the_header(self):
        parsed = measurement_points(["Layers", "1", "2", "3"], ["Voltage", "456", "418", "384"])
        self.assertEqual(parsed["points"], [(1.0, 456.0), (2.0, 418.0), (3.0, 384.0)])
        self.assertEqual(parsed["skipped_rows"], 1)
        self.assertIsNone(parsed["no_film_voltage_mv"])

    def test_points_are_sorted_by_thickness(self):
        parsed = measurement_points(["3", "1", "2"], ["384", "456", "418"])
        self.assertEqual([value for value, _ in parsed["points"]], [1.0, 2.0, 3.0])

    def test_units_and_thousands_separators_are_tolerated(self):
        parsed = measurement_points(["1 layer", "2 layers"], ["1,456 mV", "1,418mV"])
        self.assertEqual(parsed["points"], [(1.0, 1456.0), (2.0, 1418.0)])

    def test_zero_thickness_is_taken_as_the_no_film_reading(self):
        parsed = measurement_points(["0", "1", "2"], ["498", "456", "418"])
        self.assertEqual(parsed["no_film_voltage_mv"], 498.0)
        self.assertEqual(parsed["points"], [(1.0, 456.0), (2.0, 418.0)])

    def test_misaligned_selections_are_refused_not_guessed(self):
        with self.assertRaises(ValueError) as caught:
            measurement_points(["1", "2", "3"], ["456", "", "384"])
        self.assertIn("line up", str(caught.exception))

    def test_selections_of_different_sizes_are_refused(self):
        with self.assertRaises(ValueError) as caught:
            measurement_points(["1", "2"], ["456"])
        self.assertIn("same number of cells", str(caught.exception))

    def test_repeated_thickness_is_refused(self):
        with self.assertRaises(ValueError) as caught:
            measurement_points(["1", "1"], ["456", "450"])
        self.assertIn("repeats", str(caught.exception))

    def test_negative_values_and_empty_selections_are_refused(self):
        for thickness, voltage in ([["-1"], ["456"]], [["1"], ["-5"]]):
            with self.assertRaises(ValueError):
                measurement_points(thickness, voltage)
        with self.assertRaises(ValueError):
            measurement_points([], [])
        with self.assertRaises(ValueError) as caught:
            measurement_points(["Layers"], ["Voltage"])
        self.assertIn("No numeric measurement pairs", str(caught.exception))

    def test_conflicting_zero_thickness_readings_are_refused(self):
        with self.assertRaises(ValueError):
            measurement_points(["0", "0", "1"], ["498", "480", "456"])


if __name__ == "__main__":
    unittest.main()
