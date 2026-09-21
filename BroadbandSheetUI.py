"""Spreadsheet viewer for picking measured data points out of a file.

Opens an .xlsx/.csv, draws the sheet as a grid, and lets the user drag over the
cells holding the thicknesses and the cells holding the measured voltages.  It
converts nothing itself: ``MeasurementSheets.measurement_points`` turns the two
selections into data points and explains anything it refuses.
"""

from __future__ import annotations

import os
import tkinter as tk
from tkinter import messagebox

import customtkinter as ctk

from MeasurementSheets import cells_in_range, column_label, measurement_points, range_label, read_tabular_file

CELL_WIDTH = 108
CELL_HEIGHT = 24
GUTTER_WIDTH = 52
GRID_BG = "#1f2430"
HEADER_BG = "#2f3646"
LINE = "#3f4757"
TEXT = "#e2e8f0"
THICKNESS_FILL = "#1d4ed8"
VOLTAGE_FILL = "#047857"
SELECT_OUTLINE = "#fbbf24"

ROLE_THICKNESS = "thickness"
ROLE_VOLTAGE = "voltage"


class SheetRangeDialog(ctk.CTkToplevel):
    """Modal picker.  ``result`` is None when cancelled, otherwise a mapping with
    ``points``, ``no_film_voltage_mv``, ``source`` and ``summary``."""

    def __init__(self, parent, path, thickness_caption="film layers / thickness"):
        super().__init__(parent)
        self.title(f"Select measured data - {os.path.basename(path)}")
        self.geometry("1040x680")
        self.minsize(820, 520)
        self.result = None
        self.thickness_caption = thickness_caption
        self.book = read_tabular_file(path)
        self.sheet_name = next(iter(self.book["sheets"]))
        self.ranges = {ROLE_THICKNESS: None, ROLE_VOLTAGE: None}
        self._anchor = None
        self._dragging = None
        self._build()
        self._draw_sheet()
        self._update_preview()
        self.transient(parent)
        self.after(100, self._grab)

    def _grab(self):
        if self.winfo_exists():
            try:
                self.grab_set()
            except tk.TclError:
                pass

    # ------------------------------------------------------------------ layout
    def _build(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        top = ctk.CTkFrame(self, fg_color="transparent")
        top.grid(row=0, column=0, sticky="ew", padx=10, pady=(10, 4))
        ctk.CTkLabel(top, text=self.book["name"], font=("Arial", 14, "bold")).pack(side="left")
        if len(self.book["sheets"]) > 1:
            self.sheet_var = ctk.StringVar(value=self.sheet_name)
            ctk.CTkOptionMenu(top, variable=self.sheet_var, values=list(self.book["sheets"]),
                              command=self._sheet_changed, width=200).pack(side="left", padx=12)
        else:
            self.sheet_var = ctk.StringVar(value=self.sheet_name)
        ctk.CTkLabel(top, text="Drag over the cells, then press the button for what they hold.",
                     text_color="#94a3b8").pack(side="left", padx=12)

        frame = ctk.CTkFrame(self)
        frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=4)
        frame.grid_columnconfigure(0, weight=1)
        frame.grid_rowconfigure(0, weight=1)
        self.canvas = tk.Canvas(frame, background=GRID_BG, highlightthickness=0)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        y_scroll = ctk.CTkScrollbar(frame, command=self.canvas.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll = ctk.CTkScrollbar(frame, orientation="horizontal", command=self.canvas.xview)
        x_scroll.grid(row=1, column=0, sticky="ew")
        self.canvas.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.canvas.bind("<Button-1>", self._on_press)
        self.canvas.bind("<B1-Motion>", self._on_drag)
        self.canvas.bind("<ButtonRelease-1>", self._on_release)
        self.canvas.bind("<MouseWheel>", lambda e: self.canvas.yview_scroll(-int(e.delta / 120), "units"))

        buttons = ctk.CTkFrame(self, fg_color="transparent")
        buttons.grid(row=2, column=0, sticky="ew", padx=10, pady=(4, 0))
        self.thickness_button = ctk.CTkButton(buttons, text=f"Selection = {self.thickness_caption}",
                                              fg_color=THICKNESS_FILL, hover_color="#1e40af", height=32,
                                              command=lambda: self._assign(ROLE_THICKNESS))
        self.thickness_button.pack(side="left", padx=(0, 8))
        self.voltage_button = ctk.CTkButton(buttons, text="Selection = measured voltage (mV)",
                                            fg_color=VOLTAGE_FILL, hover_color="#059669", height=32,
                                            command=lambda: self._assign(ROLE_VOLTAGE))
        self.voltage_button.pack(side="left")
        ctk.CTkButton(buttons, text="Clear", width=80, height=32, fg_color="#555555",
                      command=self._clear).pack(side="left", padx=8)

        self.preview = ctk.CTkLabel(self, text="", justify="left", anchor="w", font=("Consolas", 12),
                                    wraplength=1000)
        self.preview.grid(row=3, column=0, sticky="ew", padx=12, pady=(8, 0))

        actions = ctk.CTkFrame(self, fg_color="transparent")
        actions.grid(row=4, column=0, sticky="ew", padx=10, pady=10)
        self.use_button = ctk.CTkButton(actions, text="Use this data", height=36, fg_color="#00aa00",
                                        hover_color="#008800", font=("Arial", 13, "bold"), command=self._confirm)
        self.use_button.pack(side="right")
        ctk.CTkButton(actions, text="Cancel", height=36, width=110, fg_color="#555555",
                      command=self._cancel).pack(side="right", padx=8)

    # -------------------------------------------------------------------- grid
    @property
    def grid_rows(self):
        return self.book["sheets"][self.sheet_name]

    def _column_count(self):
        return max((len(row) for row in self.grid_rows), default=1) or 1

    def _sheet_changed(self, name):
        self.sheet_name = name
        self.ranges = {ROLE_THICKNESS: None, ROLE_VOLTAGE: None}
        self._draw_sheet()
        self._update_preview()

    def _draw_sheet(self):
        canvas = self.canvas
        canvas.delete("all")
        rows, columns = len(self.grid_rows), self._column_count()
        width = GUTTER_WIDTH + columns * CELL_WIDTH
        height = CELL_HEIGHT + rows * CELL_HEIGHT
        canvas.configure(scrollregion=(0, 0, width, height))

        canvas.create_rectangle(0, 0, width, CELL_HEIGHT, fill=HEADER_BG, outline=LINE)
        canvas.create_rectangle(0, 0, GUTTER_WIDTH, height, fill=HEADER_BG, outline=LINE)
        for column in range(columns):
            x = GUTTER_WIDTH + column * CELL_WIDTH
            canvas.create_line(x, 0, x, height, fill=LINE)
            canvas.create_text(x + CELL_WIDTH / 2, CELL_HEIGHT / 2, text=column_label(column),
                               fill=TEXT, font=("Arial", 10, "bold"))
        for row_index, row in enumerate(self.grid_rows):
            y = CELL_HEIGHT + row_index * CELL_HEIGHT
            canvas.create_line(0, y, width, y, fill=LINE)
            canvas.create_text(GUTTER_WIDTH / 2, y + CELL_HEIGHT / 2, text=str(row_index + 1),
                               fill="#94a3b8", font=("Arial", 10))
            for column_index in range(columns):
                text = row[column_index] if column_index < len(row) else ""
                if not text:
                    continue
                canvas.create_text(GUTTER_WIDTH + column_index * CELL_WIDTH + 6, y + CELL_HEIGHT / 2,
                                   text=text[:16], anchor="w", fill=TEXT, font=("Consolas", 11), tags="cell")
        self._redraw_highlights()

    def _cell_at(self, event):
        x = self.canvas.canvasx(event.x) - GUTTER_WIDTH
        y = self.canvas.canvasy(event.y) - CELL_HEIGHT
        if x < 0:
            return None
        column = int(x // CELL_WIDTH)
        # Clicking a column header selects that whole column, which is the common case.
        row = -1 if y < 0 else int(y // CELL_HEIGHT)
        if column >= self._column_count() or row >= len(self.grid_rows):
            return None
        return row, column

    def _on_press(self, event):
        cell = self._cell_at(event)
        if cell is None:
            return
        row, column = cell
        if row < 0:
            self._anchor = (0, column)
            self._dragging = (len(self.grid_rows) - 1, column)
        else:
            self._anchor = self._dragging = (row, column)
        self._redraw_highlights()

    def _on_drag(self, event):
        if self._anchor is None:
            return
        cell = self._cell_at(event)
        if cell is not None and cell[0] >= 0:
            self._dragging = cell
            self._redraw_highlights()

    def _on_release(self, _event):
        self._redraw_highlights()

    def _selection(self):
        if self._anchor is None or self._dragging is None:
            return None
        (r0, c0), (r1, c1) = self._anchor, self._dragging
        return min(r0, r1), min(c0, c1), max(r0, r1), max(c0, c1)

    def _rect(self, top, left, bottom, right):
        return (GUTTER_WIDTH + left * CELL_WIDTH, CELL_HEIGHT + top * CELL_HEIGHT,
                GUTTER_WIDTH + (right + 1) * CELL_WIDTH, CELL_HEIGHT + (bottom + 1) * CELL_HEIGHT)

    def _redraw_highlights(self):
        self.canvas.delete("highlight")
        for role, fill in ((ROLE_THICKNESS, THICKNESS_FILL), (ROLE_VOLTAGE, VOLTAGE_FILL)):
            area = self.ranges[role]
            if area:
                self.canvas.create_rectangle(*self._rect(*area), fill=fill, outline=fill, stipple="gray50",
                                             tags="highlight")
        area = self._selection()
        if area:
            self.canvas.create_rectangle(*self._rect(*area), outline=SELECT_OUTLINE, width=2, tags="highlight")
        self.canvas.tag_lower("highlight", "cell")

    # ------------------------------------------------------------- assignment
    def _assign(self, role):
        area = self._selection()
        if area is None:
            messagebox.showinfo("Nothing selected", "Drag over the cells first, then press this button.", parent=self)
            return
        self.ranges[role] = area
        self._redraw_highlights()
        self._update_preview()

    def _clear(self):
        self.ranges = {ROLE_THICKNESS: None, ROLE_VOLTAGE: None}
        self._anchor = self._dragging = None
        self._redraw_highlights()
        self._update_preview()

    def parsed(self):
        """(points_mapping, error_text).  Exactly one of the two is None."""
        if not self.ranges[ROLE_THICKNESS] or not self.ranges[ROLE_VOLTAGE]:
            return None, "Select the cells for both columns."
        try:
            thickness = cells_in_range(self.grid_rows, *self.ranges[ROLE_THICKNESS])
            voltage = cells_in_range(self.grid_rows, *self.ranges[ROLE_VOLTAGE])
            return measurement_points(thickness, voltage, self.thickness_caption, "voltage"), None
        except ValueError as exc:
            return None, str(exc)

    def _range_text(self, role):
        area = self.ranges[role]
        return "—" if area is None else range_label(*area)

    def _update_preview(self):
        parsed, error = self.parsed()
        header = (f"{self.thickness_caption}: {self._range_text(ROLE_THICKNESS)}     "
                  f"voltage: {self._range_text(ROLE_VOLTAGE)}")
        if parsed is None:
            self.preview.configure(text=f"{header}\n{error}", text_color="#fca5a5")
            self.use_button.configure(state="disabled")
            return
        pairs = "   ".join(f"{value:g} -> {voltage:g} mV" for value, voltage in parsed["points"][:8])
        extra = f"   (+{len(parsed['points']) - 8} more)" if len(parsed["points"]) > 8 else ""
        lines = [header, f"{len(parsed['points'])} point(s):  {pairs}{extra}"]
        if parsed["no_film_voltage_mv"] is not None:
            lines.append(f"Zero thickness reading {parsed['no_film_voltage_mv']:g} mV will be used as the "
                         f"no-film voltage V0.")
        if parsed["skipped_rows"]:
            lines.append(f"{parsed['skipped_rows']} header/blank row(s) skipped.")
        self.preview.configure(text="\n".join(lines), text_color="#a7f3d0")
        self.use_button.configure(state="normal")

    # ---------------------------------------------------------------- closing
    def _confirm(self):
        parsed, error = self.parsed()
        if parsed is None:
            messagebox.showerror("Cannot use this selection", error, parent=self)
            return
        source = (f"{self.book['name']} [{self.sheet_name}] "
                  f"{self._range_text(ROLE_THICKNESS)} / {self._range_text(ROLE_VOLTAGE)}")
        self.result = dict(parsed, source=source,
                           summary=f"{len(parsed['points'])} point(s) from {source}")
        self.grab_release()
        self.destroy()

    def _cancel(self):
        self.result = None
        self.grab_release()
        self.destroy()


def choose_measured_points(parent, path, thickness_caption="film layers / thickness"):
    """Open the picker modally and return its result (None when cancelled)."""
    dialog = SheetRangeDialog(parent, path, thickness_caption)
    parent.wait_window(dialog)
    return dialog.result
