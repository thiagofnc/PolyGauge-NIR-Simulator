import customtkinter as ctk
import numpy as np
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import os
import itertools
import pandas as pd
from scipy.optimize import nnls

# Import the backend components
from Components import LightSource, OpticalFilter, Sensor, MaterialLayer
from Simulation import run_simulation
from ComponentDatabase import load_component_database
from ChannelOptimizer import material_names_from_stack, rank_orthogonal_combinations
from tkinter import filedialog
from MeasuredData import (LOG_DIR, discover_measured_samples, external_group_key,
                          infer_mode, load_log_spectrum, strip_mode_tokens,
                          wavenumber_to_nm)

# Categorical series palette, stepped for the dark chart surface (#2b2b2b) and
# assigned to samples in this fixed order. Past eight samples the hues repeat
# with a different line style, so every series stays a unique colour+dash pair.
MEASURED_PALETTE = ["#3987e5", "#d95926", "#199e70", "#c98500",
                    "#d55181", "#22c55e", "#9085e9", "#e66767"]
MEASURED_LINESTYLES = ["-", "--", ":"]
MEASURED_LINESTYLE_GLYPHS = {"-": "solid", "--": "dashed", ":": "dotted"}

def load_csv_spectrum(filepath, master_wl, x_type='nm', y_type='transmission', nominal_thickness_mm=0.05):
    """
    Loads raw pixel data from CSV, converts units to nm and alpha, 
    and maps it to the simulation's master wavelength array.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))
    absolute_filepath = os.path.join(base_dir, filepath)

    if not os.path.exists(absolute_filepath):
        print(f"Warning: Could not find '{absolute_filepath}'.")
        return np.zeros_like(master_wl)

    try:
        # Use the absolute_filepath here instead of filepath
        df = pd.read_csv(absolute_filepath, skipinitialspace=True)
        # Check for column names, otherwise fall back to indexes
        if 'x' in df.columns and 'y' in df.columns:
            x = df['x'].values
            y = df['y'].values
        else:
            x = df.iloc[:, 0].values
            y = df.iloc[:, 1].values
    except Exception as e:
        print(f"Error loading {absolute_filepath}: {e}")
        return np.zeros_like(master_wl)

    # 1. Convert X-axis to Nanometers
    if x_type == 'cm^-1':
        wl_nm = 1e7 / x
    else: # 'nm'
        wl_nm = x

    # Sort arrays because numpy interpolation requires strictly increasing X values
    sort_idx = np.argsort(wl_nm)
    wl_nm = wl_nm[sort_idx]
    y = y[sort_idx]

    # 2. Convert Y-axis to Absorption Coefficient (alpha) in mm^-1
    if y_type == 'transmission':
        # Beer-Lambert: Transmission T = exp(-alpha * thickness) -> alpha = -ln(T) / thickness
        y_clipped = np.clip(y, 1e-6, 1.0) # Prevent math domain errors from 0 or negative values
        alpha = -np.log(y_clipped) / nominal_thickness_mm
    else:
        # If it's already an Intensity/Absorbance measure, use it directly.
        alpha = y

    # 3. Interpolate onto the master simulation wavelength array
    # We use left=0.0, right=0.0 so we don't extrapolate garbage outside the dataset's bounds
    interp_alpha = np.interp(master_wl, wl_nm, alpha, left=0.0, right=0.0)
    return interp_alpha

# --- Physics & Data Helpers ---
def gaussian_peak(wl, center, width, amplitude):
    return amplitude * np.exp(-((wl - center) ** 2) / (2 * width ** 2))

def gaussian_bandpass(wl, center, fwhm, peak=1.0):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return gaussian_peak(wl, center, sigma, peak)

def band_model(wl, bands):
    alpha = np.zeros_like(wl, dtype=float)
    for center, width, amplitude in bands:
        alpha += gaussian_peak(wl, center, width, amplitude)
    return alpha

def blackbody_spectrum(wl_nm, temp_k):
    h, c, k = 6.626e-34, 3.0e8, 1.38e-23
    wl_m = wl_nm * 1e-9
    intensity = (2 * h * c**2 / wl_m**5) / (np.exp(h * c / (wl_m * k * temp_k)) - 1)
    return intensity / np.max(intensity)

def bounded_blackbody_spectrum(wl_nm, temp_k, min_nm=None, max_nm=None):
    spectrum = blackbody_spectrum(wl_nm, temp_k)
    mask = np.ones_like(wl_nm, dtype=bool)
    if min_nm is not None:
        mask &= wl_nm >= min_nm
    if max_nm is not None:
        mask &= wl_nm <= max_nm
    return np.where(mask, spectrum, 0.0)

def ingaas_responsivity(wl):
    return np.maximum(np.interp(wl, [1500, 2000, 2300, 2500, 2600], [0.6, 0.9, 1.0, 0.8, 0.0]), 0)

def mct_responsivity(wl):
    return np.maximum(np.interp(wl, [2500, 2800, 3200, 3800, 4100], [0.0, 0.8, 1.0, 0.9, 0.0]), 0)

def inassb_responsivity(wl):
    return np.maximum(np.interp(wl, [2600, 2700, 3300, 5000, 5300, 5400], [0.0, 0.75, 1.0, 0.95, 0.5, 0.0]), 0)

def alpha_from_k(wl_nm, k):
    wl_mm = wl_nm * 1e-6
    with np.errstate(divide='ignore', invalid='ignore'):
        return (4 * np.pi * k) / wl_mm

def load_database_file(filepath, master_wl):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    absolute_filepath = os.path.join(base_dir, filepath)

    if not os.path.exists(absolute_filepath):
        print(f"Warning: Could not find '{absolute_filepath}'.")
        return None, None

    raw_wl, raw_n, raw_k = [], [], []

    with open(absolute_filepath, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    wl = float(parts[0])
                    n = float(parts[1])
                    raw_wl.append(wl)
                    raw_n.append(n)
                    if len(parts) >= 3:
                        raw_k.append(float(parts[2]))
                except ValueError:
                    pass

    if not raw_wl:
        return None, None

    raw_wl = np.array(raw_wl) * 1000.0
    raw_n = np.array(raw_n)
    interp_n = np.interp(master_wl, raw_wl, raw_n, left=raw_n[0], right=raw_n[-1])

    if raw_k:
        raw_k = np.array(raw_k)
        interp_k = np.interp(master_wl, raw_wl, raw_k, left=raw_k[0], right=raw_k[-1])
    else:
        interp_k = None

    return interp_n, interp_k


MATERIAL_SOURCE_NOTES = {
    "PE": {
        "accuracy": "High from 2.0-4.0 um where PE_data.yml supplies n,k for LLDPE; moderate from 1.0-2.0 um where alpha is a literature-assignment band model.",
        "sources": [
            "David et al. 2022 / refractiveindex.info PE David dataset, n,k 2.00-12.3 um.",
            "Mizushima et al. 2012, PE NIR assignments in the 1650-1900 nm region.",
            "Shi et al. 2012, CH2-rich polyethylene peaks near 1210, 1730, and 1760 nm.",
            "General NIR assignment tables: aliphatic C-H first overtone 1700-1800 nm and combination 2300-2400 nm.",
        ],
    },
    "EVOH": {
        "accuracy": "Moderate for peak positions; low-to-moderate for absolute alpha because no open tabulated EVOH n,k/alpha curve was found for 1.0-4.0 um.",
        "sources": [
            "Iwamoto et al. 2001, EVOH FT-NIR OH combination at about 4780 cm^-1 and shoulder/peak at 4970 cm^-1.",
            "Su et al. 2015 and EVOH FTIR studies for hydrogen-bonded OH and CH stretching fundamentals.",
            "General NIR assignment tables for O-H overtone/combination and aliphatic C-H overtones/combinations.",
        ],
    },
    "Nylon 6": {
        "accuracy": "Moderate for peak positions; low-to-moderate for absolute alpha because open tabulated Nylon 6 n,k/alpha data was not found for this exact range.",
        "sources": [
            "Noda et al. 2016, Nylon 6 NIR bands at 1485 nm and 1535 nm for amorphous/crystalline structure.",
            "Noda et al. 2018, PA6 CH2 combination bands at 2300 nm and 2355 nm.",
            "Wu and Siesler 1999, polyamide NIR overtone/combination assignments transferable to aliphatic polyamides.",
            "General IR/NIR assignment tables for N-H, O-H, and C-H overtone/combination regions.",
        ],
    },
    "Water": {
        "accuracy": "High for pure liquid water where Water_data.yml supplies Hale-Querry n,k at 25 C.",
        "sources": [
            "Hale and Querry 1973 / refractiveindex.info water dataset, n,k 0.2-200 um at 25 C.",
        ],
    },
}

# Set UI Theme
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

class WebGaugingApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("PolyGauge NIR Simulator")
        self.geometry("1400x850")

        self.wl = np.arange(1000, 4000, 1.0)
        self.web_layers = []
        self.sensor_channels = []
        self.component_database = load_component_database()
        self.colors = ['#00ffcc', '#ff3366', '#ffcc00', '#cc33ff', '#33ccff']

        self.material_library = {
            "Air": {"alpha": np.zeros_like(self.wl), "n": 1.0, "accuracy": "Exact for a non-absorbing reference medium in this simplified model.", "sources": []},
            "PE": {
                "alpha": band_model(self.wl, [
                    (1210, 34, 0.16),  # CH2 second overtone.
                    (1410, 58, 0.12),  # CH2 stretch + deformation combination envelope.
                    (1690, 24, 0.18),  # CH3 branch/end-group contribution in LDPE/LLDPE.
                    (1710, 24, 0.22),  # CH3 chain-end contribution in LDPE/LLDPE.
                    (1730, 23, 1.00),  # CH2 first-overtone dominant band.
                    (1760, 25, 0.65),  # CH2 first-overtone shoulder.
                    (2310, 32, 2.20),  # Aliphatic C-H stretch + bend combination.
                    (2350, 34, 1.40),  # Aliphatic C-H stretch + bend combination shoulder.
                ]),
                "n": 1.51,
                **MATERIAL_SOURCE_NOTES["PE"],
            },
            "EVOH": {
                "alpha": band_model(self.wl, [
                    (1410, 42, 0.55),  # O-H first overtone / hydrogen-bond-sensitive envelope.
                    (1730, 28, 0.35),  # Ethylene CH2 first overtone.
                    (1760, 30, 0.25),  # Ethylene CH2 first-overtone shoulder.
                    (2012, 42, 0.90),  # 4970 cm^-1 OH combination shoulder/peak.
                    (2092, 55, 1.90),  # 4780 cm^-1 aggregated-OH combination.
                    (2310, 36, 0.55),  # Ethylene CH stretch + bend combination.
                    (2350, 38, 0.38),  # Ethylene CH combination shoulder.
                    (3003, 105, 18.0), # 3330 cm^-1 hydrogen-bonded O-H stretch fundamental.
                    (3409, 80, 9.5),   # 2933 cm^-1 CH2 asymmetric stretch fundamental.
                    (3506, 82, 7.0),   # 2852 cm^-1 CH2 symmetric stretch fundamental.
                ]),
                "n": 1.53,
                **MATERIAL_SOURCE_NOTES["EVOH"],
            },
            "Nylon 6": {
                "alpha": band_model(self.wl, [
                    (1485, 28, 0.85),  # N-H first overtone, amorphous PA6 contribution.
                    (1535, 30, 1.20),  # N-H first overtone, crystalline PA6 contribution.
                    (1940, 70, 0.18),  # Residual moisture sensitivity; intentionally weak for dry PA6.
                    (2040, 44, 1.55),  # N-H stretch + amide/combination region.
                    (2300, 35, 0.60),  # CH2 combination from crystalline/rigid chains.
                    (2355, 35, 0.62),  # CH2 combination from amorphous/rubbery chains.
                    (3039, 76, 18.0),  # 3291 cm^-1 N-H stretching fundamental.
                    (3410, 82, 8.0),   # CH2 asymmetric stretch fundamental.
                    (3505, 86, 6.0),   # CH2 symmetric stretch fundamental.
                ]),
                "n": 1.53,
                **MATERIAL_SOURCE_NOTES["Nylon 6"],
            },
            # "Nylon 66": {
            #     "alpha": band_model(self.wl, [(1515, 32, 0.95), (1565, 34, 1.10), (2025, 42, 1.25), (2075, 45, 1.45), (2295, 35, 0.65), (2350, 36, 0.70), (3060, 78, 17.0), (3415, 82, 8.5), (3505, 88, 6.5)]),
            #     "n": 1.54,
            # },
        }

        pe_n, pe_k = load_database_file("PE_data.yml", self.wl)
        if pe_k is not None:
            measured_region = self.wl >= 2000
            pe_alpha = self.material_library["PE"]["alpha"].copy()
            pe_alpha[measured_region] = alpha_from_k(self.wl, pe_k)[measured_region]
            pe_n_combined = np.full_like(self.wl, 1.51, dtype=float)
            if pe_n is not None:
                pe_n_combined[measured_region] = pe_n[measured_region]
            self.material_library["PE"] = {"alpha": pe_alpha, "n": pe_n_combined, **MATERIAL_SOURCE_NOTES["PE"]}

        water_n, water_k = load_database_file("Water_data.yml", self.wl)
        if water_k is not None:
            self.material_library["Water"] = {"alpha": alpha_from_k(self.wl, water_k), "n": water_n if water_n is not None else 1.33, **MATERIAL_SOURCE_NOTES["Water"]}
        else:
            self.material_library["Water"] = {
                "alpha": band_model(self.wl, [(1450, 55, 12.0), (1940, 70, 40.0)]),
                "n": 1.33,
                "accuracy": "Moderate fallback only; Water_data.yml was not available.",
                "sources": MATERIAL_SOURCE_NOTES["Water"]["sources"],
            }


        # --- LOAD CUSTOM EXTRACTED CSV DATA ---
        
        pe_csv = load_csv_spectrum("PETransissionPercentVsWavelength_in_nm_2000nm_to_5000nm.csv", 
                                   self.wl, x_type='nm', y_type='transmission', nominal_thickness_mm=0.05)
        if np.any(pe_csv):
            self.material_library["PE (Extracted)"] = {
                "alpha": pe_csv,
                "n": 1.51,
                "accuracy": "Low-to-moderate: digitized local transmission curve with assumed 0.05 mm thickness; useful as an empirical comparison, not a certified optical-constant dataset.",
                "sources": ["Local CSV: PETransissionPercentVsWavelength_in_nm_2000nm_to_5000nm.csv"],
            }

        evoh_trans = load_csv_spectrum("EVOHTransmissionPercentvsWavelength_in_nm_2000nm_to_5000nm.csv", 
                                       self.wl, x_type='nm', y_type='transmission', nominal_thickness_mm=0.015)
        if np.any(evoh_trans):
            self.material_library["EVOH 2-5um (Extracted)"] = {
                "alpha": evoh_trans,
                "n": 1.53,
                "accuracy": "Low-to-moderate: digitized local transmission curve with assumed 0.015 mm thickness; strongest peaks are plausible but absolute alpha depends on the true film thickness and digitization.",
                "sources": ["Local CSV: EVOHTransmissionPercentvsWavelength_in_nm_2000nm_to_5000nm.csv"],
            }

        evoh_int_5k = load_csv_spectrum("EVOH_Intensity_vs_Wavenumber_in_inverse_cm_5555nmto9090nm.csv", 
                                        self.wl, x_type='cm^-1', y_type='intensity')
        if np.any(evoh_int_5k):
            self.material_library["EVOH 5-9um (Extracted)"] = {
                "alpha": evoh_int_5k,
                "n": 1.53,
                "accuracy": "Low: local intensity trace is not calibrated to absorption coefficient.",
                "sources": ["Local CSV: EVOH_Intensity_vs_Wavenumber_in_inverse_cm_5555nmto9090nm.csv"],
            }

        evoh_int_3k = load_csv_spectrum("EVOH_Intensity_vs_Wavenumber_in_inverse_cm_3225nmto4000nm.csv", 
                                        self.wl, x_type='cm^-1', y_type='intensity')
        if np.any(evoh_int_3k):
            self.material_library["EVOH 3-4um (Extracted)"] = {
                "alpha": evoh_int_3k,
                "n": 1.53,
                "accuracy": "Low: local intensity trace is not calibrated to absorption coefficient.",
                "sources": ["Local CSV: EVOH_Intensity_vs_Wavenumber_in_inverse_cm_3225nmto4000nm.csv"],
            }

        self.material_display_vars = {
            name: ctk.BooleanVar(value=(name != "Air"))
            for name in self.material_library
        }

        self.setup_ui()
        self.run_live_simulation()

    def setup_ui(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=3)
        self.grid_rowconfigure(0, weight=1)

        # PANEL 1: SOURCE & SENSOR STACK
        self.left_panel = ctk.CTkFrame(self)
        self.left_panel.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        ctk.CTkLabel(self.left_panel, text="Light Source", font=("Arial", 18, "bold")).pack(pady=(10, 0))
        self.source_frame = ctk.CTkFrame(self.left_panel, fg_color="#2b2b2b")
        self.source_frame.pack(fill="x", padx=10, pady=5)

        self.source_type_var = ctk.StringVar(value="Blackbody (Halogen)")
        ctk.CTkOptionMenu(self.source_frame, variable=self.source_type_var,
                          values=["Blackbody (Halogen)", "MTE6114W-WRC LED (1460nm)", "HPIR104 Thermal Emitter", "Flat Emission (Ideal)"]).pack(pady=5)

        temp_frame = ctk.CTkFrame(self.source_frame, fg_color="transparent")
        temp_frame.pack(fill="x", pady=5)
        ctk.CTkLabel(temp_frame, text="Temp (K):").pack(side="left", padx=5)
        self.source_temp_var = ctk.StringVar(value="3000")
        ctk.CTkEntry(temp_frame, textvariable=self.source_temp_var, width=60).pack(side="left", padx=5)

        ctk.CTkLabel(self.left_panel, text="Sensor Channels", font=("Arial", 18, "bold")).pack(pady=(20, 0))
        self.sensors_container = ctk.CTkScrollableFrame(self.left_panel)
        self.sensors_container.pack(fill="both", expand=True, padx=5, pady=5)

        btn_frame = ctk.CTkFrame(self.left_panel, fg_color="transparent")
        btn_frame.pack(pady=10, fill="x", padx=5)
        ctk.CTkButton(btn_frame, text="+ Add Sensor Channel", command=self.add_sensor_ui).pack(side="left", padx=5, expand=True)
        ctk.CTkButton(btn_frame, text="Add All from Database", command=self.add_sensors_from_db, fg_color="#aa5500", hover_color="#cc7700").pack(side="left", padx=5, expand=True)

        # PANEL 2: WEB STACK
        self.stack_panel = ctk.CTkFrame(self)
        self.stack_panel.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")

        ctk.CTkLabel(self.stack_panel, text="Web Material Stack", font=("Arial", 18, "bold")).pack(pady=10)
        self.layers_container = ctk.CTkScrollableFrame(self.stack_panel)
        self.layers_container.pack(fill="both", expand=True, padx=5, pady=5)

        ctk.CTkButton(self.stack_panel, text="+ Add Layer", command=self.add_layer_ui).pack(pady=10)
        ctk.CTkButton(self.stack_panel, text="Compare Absorbance Curves", command=self.show_all_material_spectra).pack(pady=(0, 10))

        ctk.CTkLabel(self.stack_panel, text="Displayed Materials", font=("Arial", 14, "bold")).pack(pady=(5, 0))
        self.material_display_container = ctk.CTkScrollableFrame(self.stack_panel, height=130)
        self.material_display_container.pack(fill="x", padx=5, pady=(5, 10))
        for mat_name, display_var in self.material_display_vars.items():
            if mat_name == "Air":
                continue
            ctk.CTkCheckBox(
                self.material_display_container,
                text=mat_name,
                variable=display_var,
                command=self.run_live_simulation,
            ).pack(anchor="w", padx=5, pady=3)

        self.add_layer_ui(mat="PE", thick="0.05")
        self.add_layer_ui(mat="EVOH", thick="0.015")
        self.add_layer_ui(mat="Nylon 6", thick="0.02")
        self.add_layer_ui(mat="Water", thick="0.005")

        # PANEL 3: GRAPH & GLOBAL CONTROLS
        self.data_panel = ctk.CTkFrame(self)
        self.data_panel.grid(row=0, column=2, padx=10, pady=10, sticky="nsew")

        self.fig, (self.ax_top, self.ax_bot) = plt.subplots(2, 1, figsize=(8, 6), facecolor='#2b2b2b', gridspec_kw={'height_ratios': [2, 1]})
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.data_panel)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

        ctk.CTkButton(self.data_panel, text="UPDATE SIMULATION", command=self.run_live_simulation,
                      fg_color="#00aa00", hover_color="#008800", font=("Arial", 16, "bold"), height=40).pack(fill="x", padx=20, pady=10)
        ctk.CTkButton(self.data_panel, text="CHANNEL MATRIX", command=self.show_channel_matrix,
                      fg_color="#0055aa", hover_color="#0077cc", font=("Arial", 15, "bold"), height=36).pack(fill="x", padx=20, pady=(0, 10))
        ctk.CTkButton(self.data_panel, text="SOLVE THICKNESSES", command=self.show_thickness_solver,
                      fg_color="#0a7c66", hover_color="#0d9479", font=("Arial", 15, "bold"), height=36).pack(fill="x", padx=20, pady=(0, 10))
        ctk.CTkButton(self.data_panel, text="RANK FILTER COMBOS", command=self.show_ranked_combinations,
                      fg_color="#6b3fa0", hover_color="#7d51bd", font=("Arial", 15, "bold"), height=36).pack(fill="x", padx=20, pady=(0, 10))
        ctk.CTkButton(self.data_panel, text="MEASURED SAMPLES", command=self.show_measured_samples,
                      fg_color="#a15c00", hover_color="#c47200", font=("Arial", 15, "bold"), height=36).pack(fill="x", padx=20, pady=(0, 10))

    # --- UI Builders ---
    def add_sensors_from_db(self):
        for filter_def in self.component_database.get("filters", []):
            center = filter_def.get("center_nm")
            width = filter_def.get("fwhm_nm")
            if center is None or width is None:
                continue

            price = filter_def.get("price_usd", "")
            name = filter_def.get("name", f"{center} nm")
            notes = filter_def.get("notes")
            if notes:
                name = f"{notes} ({name})"

            self.add_sensor_ui(center=str(center), width=str(width), name=name, price=str(price))

    def add_layer_ui(self, mat="PE", thick="0.05"):
        if mat not in self.material_library:
            mat = "PE"

        frame = ctk.CTkFrame(self.layers_container, fg_color="#333333")
        frame.pack(fill="x", pady=5, padx=5)

        mat_var = ctk.StringVar(value=mat)
        ctk.CTkOptionMenu(frame, variable=mat_var, values=list(self.material_library.keys()), width=90).pack(side="left", padx=5, pady=10)

        thick_var = ctk.StringVar(value=thick)
        ctk.CTkEntry(frame, textvariable=thick_var, width=60).pack(side="left", padx=5)
        ctk.CTkLabel(frame, text="mm").pack(side="left")

        ctk.CTkButton(frame, text="X", width=30, fg_color="#aa0000", hover_color="#ff0000",
                      command=lambda f=frame: self.remove_element(f, self.web_layers)).pack(side="right", padx=(5, 5))

        ctk.CTkButton(frame, text="👁️", width=30, fg_color="#0055aa", hover_color="#0077ff",
                      command=lambda m=mat_var: self.show_material_spectra(m.get())).pack(side="right", padx=(5, 0))

        self.web_layers.append({"frame": frame, "mat_var": mat_var, "thick_var": thick_var})

    def add_sensor_ui(self, center="2310", width="15", name="Sensor", sensor_type="InGaAs", price="150"):
        color = self.colors[len(self.sensor_channels) % len(self.colors)]
        frame = ctk.CTkFrame(self.sensors_container, border_width=2, border_color=color, fg_color="#333333")
        frame.pack(fill="x", pady=5, padx=5)

        row1 = ctk.CTkFrame(frame, fg_color="transparent")
        row1.pack(fill="x", padx=5, pady=(5,0))
        name_var = ctk.StringVar(value=name)
        ctk.CTkEntry(row1, textvariable=name_var, width=180).pack(side="left")

        sensor_type_var = ctk.StringVar(value=sensor_type)
        ctk.CTkOptionMenu(row1, variable=sensor_type_var,
                          values=["InGaAs", "InAsSb (2.7-5.3um)", "MCT (MIR)", "Ideal (Flat)"],
                          width=110).pack(side="right")

        row2 = ctk.CTkFrame(frame, fg_color="transparent")
        row2.pack(fill="x", padx=5, pady=5)
        ctk.CTkLabel(row2, text="CWL:").pack(side="left")
        center_var = ctk.StringVar(value=center)
        ctk.CTkEntry(row2, textvariable=center_var, width=45).pack(side="left", padx=(0,5))

        ctk.CTkLabel(row2, text="FWHM:").pack(side="left")
        width_var = ctk.StringVar(value=width)
        ctk.CTkEntry(row2, textvariable=width_var, width=35).pack(side="left")

        ctk.CTkLabel(row2, text="Price: $").pack(side="left", padx=(5,0))
        price_var = ctk.StringVar(value=str(price))
        ctk.CTkEntry(row2, textvariable=price_var, width=40).pack(side="left")

        row3 = ctk.CTkFrame(frame, fg_color="transparent")
        row3.pack(fill="x", padx=5, pady=(0,5))
        lbl_readout = ctk.CTkLabel(row3, text="Signal: 0.00", font=("Arial", 14, "bold"), text_color=color)
        lbl_readout.pack(side="left")

        ctk.CTkButton(row3, text="X", width=30, fg_color="#aa0000", hover_color="#ff0000",
                      command=lambda f=frame: self.remove_element(f, self.sensor_channels)).pack(side="right")

        self.sensor_channels.append({
            "frame": frame, "name_var": name_var, "center_var": center_var, "width_var": width_var,
            "price_var": price_var, "sensor_type_var": sensor_type_var, "lbl_readout": lbl_readout, "color": color
        })

    def remove_element(self, frame, target_list):
        frame.destroy()
        target_list[:] = [item for item in target_list if item["frame"] != frame]

    def get_source_spectra(self):
        source_type = self.source_type_var.get()
        if source_type == "Blackbody (Halogen)":
            try: temp = float(self.source_temp_var.get())
            except ValueError: temp = 3000
            return blackbody_spectrum(self.wl, temp)
        if source_type == "MTE6114W-WRC LED (1460nm)":
            return gaussian_bandpass(self.wl, 1460, 103, 1.0)
        if source_type == "HPIR104 Thermal Emitter":
            return bounded_blackbody_spectrum(self.wl, 903.15, min_nm=2000, max_nm=11000)
        return np.ones_like(self.wl)

    def get_sensor_spectra(self, sensor_type):
        if sensor_type == "InGaAs": return ingaas_responsivity(self.wl)
        if sensor_type == "InAsSb (2.7-5.3um)": return inassb_responsivity(self.wl)
        if sensor_type == "MCT (MIR)": return mct_responsivity(self.wl)
        return np.ones_like(self.wl)

    def get_channel_definitions(self):
        channels = []
        for ch in self.sensor_channels:
            try:
                center = float(ch["center_var"].get())
                width = float(ch["width_var"].get())
                price = float(ch["price_var"].get().replace('$', '').strip())
            except ValueError:
                continue

            channels.append({
                "name": ch["name_var"].get(),
                "center": center,
                "width": width,
                "price": price,
                "sensor_type": ch["sensor_type_var"].get(),
            })

        return channels

    def build_effective_channel_matrix(self, channels, materials):
        source_spectra = self.get_source_spectra()
        matrix = np.zeros((len(channels), len(materials)), dtype=float)
        channel_weights = np.zeros(len(channels), dtype=float)

        for row_idx, channel in enumerate(channels):
            filter_spectra = gaussian_bandpass(self.wl, channel["center"], channel["width"], 1.0)
            sensor_spectra = self.get_sensor_spectra(channel["sensor_type"])
            weight = source_spectra * filter_spectra * sensor_spectra
            weight_area = np.trapezoid(weight, self.wl)
            channel_weights[row_idx] = weight_area

            if weight_area <= 0:
                continue

            for col_idx, material in enumerate(materials):
                alpha = np.asarray(self.material_library[material]["alpha"], dtype=float)
                matrix[row_idx, col_idx] = np.trapezoid(weight * alpha, self.wl) / weight_area

        return matrix, channel_weights

    def build_web_stack_objects(self):
        web_stack_objs = []
        for i, layer_data in enumerate(self.web_layers):
            mat = layer_data["mat_var"].get()
            try:
                thickness = float(layer_data["thick_var"].get())
            except ValueError:
                thickness = 0.0

            web_stack_objs.append(MaterialLayer(
                name=f"Layer_{i}_{mat}",
                thickness=thickness,
                raw_wavelengths=self.wl,
                alpha_data=self.material_library[mat]["alpha"],
                n_data=self.material_library[mat]["n"],
            ))

        return web_stack_objs

    def get_unique_stack_materials(self):
        return [
            material for material in material_names_from_stack(self.material_library, self.web_layers)
            if self.is_material_displayed(material)
        ]

    def get_displayed_material_names(self, include_air=False):
        displayed = []
        for name in self.material_library:
            if not include_air and name == "Air":
                continue

            display_var = self.material_display_vars.get(name)
            if display_var is None or display_var.get():
                displayed.append(name)

        return displayed

    def is_material_displayed(self, mat_name):
        if mat_name == "Air":
            return True
        display_var = self.material_display_vars.get(mat_name)
        return display_var is None or display_var.get()

    def show_material_spectra(self, mat_name):
        if mat_name not in self.material_library or not self.is_material_displayed(mat_name):
            return

        data = self.material_library[mat_name]
        popup = ctk.CTkToplevel(self)
        popup.title(f"Material Profile: {mat_name}")
        popup.geometry("600x450")
        popup.attributes("-topmost", True)

        fig, ax_alpha = plt.subplots(figsize=(6, 4), facecolor='#2b2b2b')
        ax_alpha.set_facecolor('#2b2b2b')
        ax_alpha.set_title(f"{mat_name} Base Spectral Properties", color='white')
        ax_alpha.set_xlabel("Wavelength (nm)", color='white')
        ax_alpha.tick_params(colors='white')

        color_alpha = '#ff3366'
        ax_alpha.set_ylabel("Absorption Coefficient (mm⁻¹)", color=color_alpha)
        ax_alpha.plot(self.wl, data["alpha"], color=color_alpha, label="Absorption (α)")
        ax_alpha.tick_params(axis='y', labelcolor=color_alpha)
        ax_alpha.grid(True, alpha=0.2)

        ax_n = ax_alpha.twinx()
        color_n = '#00ccff'
        ax_n.set_ylabel("Refractive Index (n)", color=color_n)

        n_data = np.full_like(self.wl, data["n"]) if isinstance(data["n"], (int, float)) else data["n"]
        ax_n.plot(self.wl, n_data, color=color_n, linestyle='--', label="Refractive Index (n)")
        ax_n.tick_params(axis='y', labelcolor=color_n)

        lines_1, labels_1 = ax_alpha.get_legend_handles_labels()
        lines_2, labels_2 = ax_n.get_legend_handles_labels()
        ax_alpha.legend(lines_1 + lines_2, labels_1 + labels_2, facecolor='#333333', edgecolor='white', labelcolor='white', loc='upper left')

        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=popup)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

        marker_alpha, = ax_alpha.plot([], [], marker="o", color="white", markersize=6, linestyle="None")
        marker_n, = ax_n.plot([], [], marker="o", color="white", markersize=6, linestyle="None")
        annotation = ax_alpha.annotate(
            "", xy=(0, 0), xytext=(12, 12), textcoords="offset points", color="white",
            bbox=dict(boxstyle="round,pad=0.3", fc="#333333", ec="white", alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="white"),
        )
        annotation.set_visible(False)

        def on_click(event):
            if event.xdata is None: return
            idx = int(np.argmin(np.abs(self.wl - event.xdata)))
            x_val = self.wl[idx]

            if event.inaxes == ax_n:
                y_val = n_data[idx]
                marker_n.set_data([x_val], [y_val])
                marker_alpha.set_data([], [])
                annotation.xy = (x_val, data["alpha"][idx])
                annotation.set_text(f"{x_val:.0f} nm\nn: {y_val:.4g}")
            elif event.inaxes == ax_alpha:
                y_val = data["alpha"][idx]
                marker_alpha.set_data([x_val], [y_val])
                marker_n.set_data([], [])
                annotation.xy = (x_val, y_val)
                annotation.set_text(f"{x_val:.0f} nm\nalpha: {y_val:.4g} mm^-1")
            else: return

            annotation.set_visible(True)
            canvas.draw_idle()

        canvas.mpl_connect("button_press_event", on_click)
        canvas.draw()

    def show_all_material_spectra(self):
        """Shows all material absorption curves in one popup, with optional filter overlays."""
        popup = ctk.CTkToplevel(self)
        popup.title("Material Absorbance Comparison")
        popup.geometry("900x700")
        popup.attributes("-topmost", True)

        # UI for toggling filters
        ctrl_frame = ctk.CTkFrame(popup, fg_color="transparent")
        ctrl_frame.pack(fill="x", padx=10, pady=(10, 0))

        show_filters_var = ctk.BooleanVar(value=True)

        fig, (ax_raw, ax_norm) = plt.subplots(2, 1, figsize=(9, 6), facecolor='#2b2b2b', sharex=True)

        colors = {"PE": "#00ffcc", "EVOH": "#ffcc00", "Nylon 6": "#cc33ff", "Nylon 66": "#ff66cc", "Water": "#33ccff", "Air": "#888888"}

        for ax in (ax_raw, ax_norm):
            ax.set_facecolor('#2b2b2b')
            ax.tick_params(colors='white')
            ax.grid(True, alpha=0.2)
            for spine in ax.spines.values(): spine.set_color('gray')

        plotted = []
        for mat_name in self.get_displayed_material_names():
            data = self.material_library[mat_name]

            alpha = np.asarray(data["alpha"], dtype=float)
            color = colors.get(mat_name, None)
            raw_line, = ax_raw.plot(self.wl, alpha, label=mat_name, color=color, linewidth=1.8)

            max_alpha = np.nanmax(alpha)
            if max_alpha > 0:
                norm_alpha = alpha / max_alpha
                norm_line, = ax_norm.plot(self.wl, norm_alpha, label=mat_name, color=color, linewidth=1.8)
                plotted.append({"name": mat_name, "raw": alpha, "norm": norm_alpha, "raw_line": raw_line, "norm_line": norm_line})

        ax_raw.set_title("Absorption Coefficient", color='white')
        ax_raw.set_ylabel("alpha (mm^-1)", color='white')
        ax_raw.set_yscale("symlog", linthresh=0.01)

        ax_norm.set_title("Normalized Band Shapes", color='white')
        ax_norm.set_xlabel("Wavelength (nm)", color='white')
        ax_norm.set_ylabel("Relative alpha", color='white')

        if plotted:
            for ax in (ax_raw, ax_norm):
                ax.legend(facecolor='#333333', edgecolor='white', labelcolor='white', loc='upper right')
        else:
            ax_raw.text(0.5, 0.5, "No materials selected", transform=ax_raw.transAxes, ha="center", va="center", color="white")
            ax_norm.text(0.5, 0.5, "Select materials in the main panel", transform=ax_norm.transAxes, ha="center", va="center", color="white")

        # --- Filter Overlay Logic ---
        ax_raw_twin = ax_raw.twinx()
        ax_raw_twin.set_ylabel("Filter Transmission", color='#aaaaaa')
        ax_raw_twin.tick_params(axis='y', colors='#aaaaaa')
        ax_raw_twin.set_ylim(0, 1.05)
        for spine in ax_raw_twin.spines.values(): spine.set_color('gray')

        filter_artists = []

        for ch in self.sensor_channels:
            try:
                c_wl = float(ch["center_var"].get())
                fwhm = float(ch["width_var"].get())
            except ValueError:
                continue

            color = ch["color"]
            filter_spectra = gaussian_bandpass(self.wl, c_wl, fwhm, 1.0)

            # Plot on secondary top axis
            f1 = ax_raw_twin.fill_between(self.wl, 0, filter_spectra, color=color, alpha=0.15)
            l1, = ax_raw_twin.plot(self.wl, filter_spectra, color=color, linestyle=":", alpha=0.8)

            # Plot on bottom normalized axis
            f2 = ax_norm.fill_between(self.wl, 0, filter_spectra, color=color, alpha=0.15)
            l2, = ax_norm.plot(self.wl, filter_spectra, color=color, linestyle=":", alpha=0.8)

            filter_artists.extend([f1, l1, f2, l2])

        def toggle_filters():
            is_visible = show_filters_var.get()
            for artist in filter_artists:
                artist.set_visible(is_visible)
            ax_raw_twin.get_yaxis().set_visible(is_visible)
            canvas.draw_idle()

        ctk.CTkCheckBox(ctrl_frame, text="Overlay Sensor Filters", variable=show_filters_var, command=toggle_filters).pack(side="left")

        # Hide twin y-axis by default if no filters exist
        if not filter_artists:
            ax_raw_twin.get_yaxis().set_visible(False)
            show_filters_var.set(False)

        # ----------------------------

        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=popup)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

        raw_marker, = ax_raw.plot([], [], marker="o", color="white", markersize=6, linestyle="None")
        norm_marker, = ax_norm.plot([], [], marker="o", color="white", markersize=6, linestyle="None")
        raw_annotation = ax_raw.annotate("", xy=(0, 0), xytext=(12, 12), textcoords="offset points", color="white", bbox=dict(boxstyle="round,pad=0.3", fc="#333333", ec="white", alpha=0.95), arrowprops=dict(arrowstyle="->", color="white"))
        norm_annotation = ax_norm.annotate("", xy=(0, 0), xytext=(12, 12), textcoords="offset points", color="white", bbox=dict(boxstyle="round,pad=0.3", fc="#333333", ec="white", alpha=0.95), arrowprops=dict(arrowstyle="->", color="white"))
        raw_annotation.set_visible(False)
        norm_annotation.set_visible(False)

        def on_click(event):
            # Treat clicks on the twin filter axis as clicks on the main raw axis
            click_ax = ax_raw if event.inaxes == ax_raw_twin else event.inaxes
            if event.xdata is None or click_ax not in (ax_raw, ax_norm): return

            idx = int(np.argmin(np.abs(self.wl - event.xdata)))
            x_val = self.wl[idx]
            series_key = "raw" if click_ax == ax_raw else "norm"
            click_xy = np.array([event.x, event.y])

            nearest = None
            nearest_distance = float("inf")
            for item in plotted:
                y_val = item[series_key][idx]
                pixel_xy = click_ax.transData.transform((x_val, y_val))
                distance = np.linalg.norm(pixel_xy - click_xy)
                if distance < nearest_distance:
                    nearest = item
                    nearest_distance = distance

            if nearest is None: return

            y_val = nearest[series_key][idx]
            if click_ax == ax_raw:
                raw_marker.set_data([x_val], [y_val])
                norm_marker.set_data([], [])
                raw_annotation.xy = (x_val, y_val)
                raw_annotation.set_text(f"{nearest['name']}\n{x_val:.0f} nm\nalpha: {y_val:.4g} mm^-1")
                raw_annotation.set_visible(True)
                norm_annotation.set_visible(False)
            else:
                norm_marker.set_data([x_val], [y_val])
                raw_marker.set_data([], [])
                norm_annotation.xy = (x_val, y_val)
                norm_annotation.set_text(f"{nearest['name']}\n{x_val:.0f} nm\nrelative: {y_val:.4g}")
                norm_annotation.set_visible(True)
                raw_annotation.set_visible(False)

            canvas.draw_idle()

        canvas.mpl_connect("button_press_event", on_click)
        canvas.draw()

    # --- Measured spectrometer data ---
    def get_measured_spectrum(self, path):
        """Loads a log file once and caches it for the session."""
        if not hasattr(self, "_measured_spectrum_cache"):
            self._measured_spectrum_cache = {}
        if path not in self._measured_spectrum_cache:
            self._measured_spectrum_cache[path] = load_log_spectrum(path)
        return self._measured_spectrum_cache[path]

    def show_measured_samples(self):
        """Plots the measured FTIR logs with a per-sample show/hide panel."""
        samples, reference_headers = discover_measured_samples()

        popup = ctk.CTkToplevel(self)
        popup.title("Measured Sample Spectra")
        popup.geometry("1350x820")
        popup.attributes("-topmost", True)

        if not samples:
            ctk.CTkLabel(
                popup,
                text="No measured logs found in logs_full_range/.",
                font=("Arial", 14, "bold"),
            ).pack(expand=True, padx=20, pady=20)
            return

        for index, sample in enumerate(samples):
            sample["color"] = MEASURED_PALETTE[index % len(MEASURED_PALETTE)]
            sample["linestyle"] = MEASURED_LINESTYLES[(index // len(MEASURED_PALETTE)) % len(MEASURED_LINESTYLES)]
            # Start with a readable handful rather than 18 overlapping traces.
            sample["visible_var"] = ctk.BooleanVar(value=index < 3)

        y_mode_var = ctk.StringVar(value="Absorbance")
        x_mode_var = ctk.StringVar(value="Wavenumber (1/cm)")
        x_min_var = ctk.StringVar(value="")
        x_max_var = ctk.StringVar(value="")
        range_status_var = ctk.StringVar(value="")
        bandwidth_var = ctk.StringVar(value="400")
        filter_status_var = ctk.StringVar(value="Click a point on the graph, then Add Filter")

        # Filter bands are stored as (low, high) wavenumber edges so they stay
        # put when the X axis switches between cm^-1 and nm.
        filter_bands = []
        selection = {"wavenumber": None, "display_x": None}

        popup.grid_columnconfigure(0, weight=1)
        popup.grid_columnconfigure(1, weight=0)
        popup.grid_rowconfigure(0, weight=1)

        plot_frame = ctk.CTkFrame(popup, fg_color="transparent")
        plot_frame.grid(row=0, column=0, sticky="nsew", padx=(10, 5), pady=10)

        fig, ax = plt.subplots(figsize=(9, 6.5), facecolor='#2b2b2b')
        ax.set_facecolor('#2b2b2b')
        ax.tick_params(colors='white')
        ax.grid(True, alpha=0.2)
        for spine in ax.spines.values():
            spine.set_color('gray')

        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)

        annotation = ax.annotate(
            "", xy=(0, 0), xytext=(12, 12), textcoords="offset points", color="white",
            bbox=dict(boxstyle="round,pad=0.3", fc="#333333", ec="white", alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="white"),
        )
        annotation.set_visible(False)
        marker, = ax.plot([], [], marker="o", color="white", markersize=6, linestyle="None")

        plotted = []
        home_limits = {}

        status_widget = {}
        filter_status_widget = {}

        def set_status(text, error=False):
            range_status_var.set(text)
            label = status_widget.get("label")
            if label is not None:
                label.configure(text_color="#e66767" if error else "#aaaaaa")

        def set_filter_status(text, error=False):
            filter_status_var.set(text)
            label = filter_status_widget.get("label")
            if label is not None:
                label.configure(text_color="#e66767" if error else "#aaaaaa")

        def get_range_bounds():
            """Parses the X-range entries into (low, high) in the current unit.

            Either side may be blank for 'unbounded'; returns None on bad input
            after reporting it in the status line.
            """
            bounds = []
            for var, name in ((x_min_var, "From"), (x_max_var, "To")):
                text = var.get().strip()
                if not text:
                    bounds.append(None)
                    continue
                try:
                    bounds.append(float(text.replace(",", ".")))
                except ValueError:
                    set_status(f"{name}: '{text}' is not a number", error=True)
                    return None
            low, high = bounds
            if low is not None and high is not None and low > high:
                low, high = high, low  # Accept the bounds in either order.
            return low, high

        def redraw():
            mode_key = "absorbance" if y_mode_var.get() == "Absorbance" else "transmittance"
            use_wavenumber = x_mode_var.get().startswith("Wavenumber")

            bounds = get_range_bounds()
            if bounds is None:
                return
            range_low, range_high = bounds

            ax.clear()
            ax.set_facecolor('#2b2b2b')
            ax.tick_params(colors='white')
            ax.grid(True, alpha=0.2)
            for spine in ax.spines.values():
                spine.set_color('gray')

            plotted.clear()
            missing = []
            selected_count = 0
            for sample in samples:
                if not sample["visible_var"].get():
                    continue
                selected_count += 1
                path = sample["paths"].get(mode_key)
                if not path:
                    missing.append(sample["label"])
                    continue
                wavenumber, values, _ = self.get_measured_spectrum(path)
                if wavenumber.size == 0:
                    continue
                x_values = wavenumber if use_wavenumber else wavenumber_to_nm(wavenumber)

                # Keep only the points inside the requested X range so the Y
                # axis scales to what is actually on screen.
                if range_low is not None or range_high is not None:
                    keep = np.ones_like(x_values, dtype=bool)
                    if range_low is not None:
                        keep &= x_values >= range_low
                    if range_high is not None:
                        keep &= x_values <= range_high
                    x_values, values = x_values[keep], values[keep]
                    if x_values.size == 0:
                        continue

                ax.plot(x_values, values, label=sample["label"], color=sample["color"],
                        linestyle=sample["linestyle"], linewidth=1.6)
                plotted.append({"name": sample["label"], "x": x_values, "y": values})

            # Filter bands sit under the traces as a recessive highlight, so
            # they never compete with the sample colours.
            for band in filter_bands:
                low_cm, high_cm = band["low_cm"], band["high_cm"]
                if use_wavenumber:
                    left, right = low_cm, high_cm
                    centre = band["centre_cm"]
                else:
                    left, right = wavenumber_to_nm(high_cm), wavenumber_to_nm(low_cm)
                    centre = wavenumber_to_nm(band["centre_cm"])
                ax.axvspan(left, right, color="#ffffff", alpha=0.10, zorder=0)
                ax.axvline(centre, color="#dddddd", alpha=0.45, linewidth=1.0,
                           linestyle="--", zorder=0)

            unit = "cm^-1" if use_wavenumber else "nm"
            if range_low is None and range_high is None:
                set_status(f"Full range ({unit})")
            else:
                low_text = f"{range_low:g}" if range_low is not None else "min"
                high_text = f"{range_high:g}" if range_high is not None else "max"
                points = plotted[0]["x"].size if plotted else 0
                set_status(f"{low_text} - {high_text} {unit}  ({points} pts/trace)",
                           error=bool(selected_count and not plotted and not missing))

            if use_wavenumber:
                ax.set_xlabel("Wavenumber (cm^-1)", color='white')
                ax.invert_xaxis()  # Conventional FTIR orientation.
            else:
                ax.set_xlabel("Wavelength (nm)", color='white')

            if mode_key == "absorbance":
                ax.set_ylabel("Absorbance (a.u.)", color='white')
                ax.set_title("Measured Absorbance", color='white')
            else:
                ax.set_ylabel("Transmittance (%)", color='white')
                ax.set_title("Measured Transmittance", color='white')

            if plotted:
                # Beyond ~10 traces the in-plot legend eats the chart; the
                # colour-coded sample list on the right serves as the legend.
                if len(plotted) <= 10:
                    ax.legend(facecolor='#333333', edgecolor='white', labelcolor='white',
                              loc='upper right', fontsize=8)
                else:
                    ax.set_title(ax.get_title() + f"  -  {len(plotted)} samples (see list at right)",
                                 color='white')
            else:
                message = "Select samples in the panel on the right"
                if missing:
                    message = "No " + y_mode_var.get().lower() + " log for the selected sample(s)"
                elif selected_count:
                    message = "No data points inside the selected X range"
                ax.text(0.5, 0.5, message, transform=ax.transAxes,
                        ha="center", va="center", color="white")

            nonlocal annotation, marker
            marker, = ax.plot([], [], marker="o", color="white", markersize=6, linestyle="None")
            annotation = ax.annotate(
                "", xy=(0, 0), xytext=(12, 12), textcoords="offset points", color="white",
                bbox=dict(boxstyle="round,pad=0.3", fc="#333333", ec="white", alpha=0.95),
                arrowprops=dict(arrowstyle="->", color="white"),
            )
            annotation.set_visible(False)

            fig.tight_layout()
            # Remember the full view so scroll-zoom can be reset to it.
            home_limits["x"] = ax.get_xlim()
            home_limits["y"] = ax.get_ylim()
            canvas.draw_idle()

        def reset_zoom():
            if home_limits.get("x"):
                ax.set_xlim(home_limits["x"])
                ax.set_ylim(home_limits["y"])
                canvas.draw_idle()

        def on_scroll(event):
            """Mouse wheel zooms about the cursor.

            Plain scroll zooms both axes; Ctrl locks the zoom to X, Shift locks
            it to Y. The arithmetic is relative to the current limits, so it
            behaves the same on the inverted wavenumber axis.
            """
            if event.inaxes != ax:
                return
            step = getattr(event, "step", 0) or (1 if event.button == "up" else -1)
            scale = 1.25 ** (-step)  # Scroll up shrinks the visible range.

            key = (event.key or "").lower()
            zoom_x = "shift" not in key
            zoom_y = "control" not in key and "ctrl" not in key

            if zoom_x and event.xdata is not None:
                left, right = ax.get_xlim()
                ax.set_xlim(event.xdata + (left - event.xdata) * scale,
                            event.xdata + (right - event.xdata) * scale)
            if zoom_y and event.ydata is not None:
                bottom, top = ax.get_ylim()
                ax.set_ylim(event.ydata + (bottom - event.ydata) * scale,
                            event.ydata + (top - event.ydata) * scale)
            canvas.draw_idle()

        canvas.mpl_connect("scroll_event", on_scroll)
        # Modifier keys only reach matplotlib when the canvas has focus.
        plot_widget = canvas.get_tk_widget()
        plot_widget.bind("<Enter>", lambda _event: plot_widget.focus_set())

        pan_state = {"active": False, "moved": False}

        def show_nearest_point(event):
            if event.xdata is None or not plotted:
                return
            click_xy = np.array([event.x, event.y])
            nearest = None
            nearest_distance = float("inf")
            for item in plotted:
                idx = int(np.argmin(np.abs(item["x"] - event.xdata)))
                pixel_xy = ax.transData.transform((item["x"][idx], item["y"][idx]))
                distance = np.linalg.norm(pixel_xy - click_xy)
                if distance < nearest_distance:
                    nearest = (item, idx)
                    nearest_distance = distance

            item, idx = nearest
            x_val, y_val = item["x"][idx], item["y"][idx]
            unit = "cm^-1" if x_mode_var.get().startswith("Wavenumber") else "nm"
            value_label = "Abs" if y_mode_var.get() == "Absorbance" else "%T"
            marker.set_data([x_val], [y_val])
            annotation.xy = (x_val, y_val)
            annotation.set_text(f"{item['name']}\n{x_val:.1f} {unit}\n{value_label}: {y_val:.4g}")
            annotation.set_visible(True)
            # Remembered in wavenumber so "Add Filter" works in either unit.
            selection["wavenumber"] = x_val if unit == "cm^-1" else 1.0e7 / x_val
            selection["display_x"] = x_val
            set_filter_status(f"Selected {x_val:.1f} {unit}")
            canvas.draw_idle()

        def on_mouse_press(event):
            if event.inaxes != ax or event.button != 1:
                return
            if event.dblclick:
                reset_zoom()
                return
            pan_state.update(
                active=True,
                moved=False,
                press_x=event.x,
                press_y=event.y,
                x_limits=ax.get_xlim(),
                y_limits=ax.get_ylim(),
            )

        def on_mouse_move(event):
            if not pan_state["active"] or event.inaxes != ax:
                return
            pixel_distance = np.hypot(event.x - pan_state["press_x"],
                                      event.y - pan_state["press_y"])
            if pixel_distance > 3:
                pan_state["moved"] = True
            left, right = pan_state["x_limits"]
            bottom, top = pan_state["y_limits"]
            x_shift = (pan_state["press_x"] - event.x) * (right - left) / ax.bbox.width
            y_shift = (pan_state["press_y"] - event.y) * (top - bottom) / ax.bbox.height
            ax.set_xlim(left + x_shift, right + x_shift)
            ax.set_ylim(bottom + y_shift, top + y_shift)
            annotation.set_visible(False)
            canvas.draw_idle()

        def on_mouse_release(event):
            if event.button != 1 or not pan_state["active"]:
                return
            pan_state["active"] = False
            if not pan_state["moved"] and event.inaxes == ax:
                show_nearest_point(event)

        canvas.mpl_connect("button_press_event", on_mouse_press)
        canvas.mpl_connect("motion_notify_event", on_mouse_move)
        canvas.mpl_connect("button_release_event", on_mouse_release)

        # --- Right-hand sample panel ---
        side_panel = ctk.CTkFrame(popup, width=380)
        side_panel.grid(row=0, column=1, sticky="nsew", padx=(5, 10), pady=10)
        side_panel.grid_propagate(False)

        ctk.CTkLabel(side_panel, text="Samples", font=("Arial", 18, "bold")).pack(pady=(10, 5))

        ctk.CTkLabel(side_panel, text="Quantity", font=("Arial", 12, "bold")).pack(pady=(5, 0))
        ctk.CTkSegmentedButton(side_panel, variable=y_mode_var,
                               values=["Absorbance", "Transmittance"],
                               command=lambda _value: redraw()).pack(fill="x", padx=10, pady=5)

        ctk.CTkLabel(side_panel, text="X axis", font=("Arial", 12, "bold")).pack(pady=(5, 0))

        x_mode_state = {"current": x_mode_var.get()}

        def on_x_mode_change(value):
            """Carries an existing range across the unit switch."""
            if value != x_mode_state["current"] and (x_min_var.get().strip() or x_max_var.get().strip()):
                bounds = get_range_bounds()
                if bounds is not None:
                    low, high = bounds
                    # nm = 1e7 / cm^-1 in both directions, and the conversion
                    # reverses the ordering, so the two ends swap.
                    new_low = f"{1.0e7 / high:.6g}" if high else ""
                    new_high = f"{1.0e7 / low:.6g}" if low else ""
                    x_min_var.set(new_low)
                    x_max_var.set(new_high)
            x_mode_state["current"] = value
            redraw()

        ctk.CTkSegmentedButton(side_panel, variable=x_mode_var,
                               values=["Wavenumber (1/cm)", "Wavelength (nm)"],
                               command=on_x_mode_change).pack(fill="x", padx=10, pady=5)

        ctk.CTkLabel(side_panel, text="X range", font=("Arial", 12, "bold")).pack(pady=(8, 0))
        range_frame = ctk.CTkFrame(side_panel, fg_color="transparent")
        range_frame.pack(fill="x", padx=10, pady=(2, 0))
        ctk.CTkLabel(range_frame, text="From", font=("Arial", 11), width=36).pack(side="left")
        min_entry = ctk.CTkEntry(range_frame, textvariable=x_min_var, width=70, placeholder_text="min")
        min_entry.pack(side="left", padx=2)
        ctk.CTkLabel(range_frame, text="To", font=("Arial", 11), width=24).pack(side="left")
        max_entry = ctk.CTkEntry(range_frame, textvariable=x_max_var, width=70, placeholder_text="max")
        max_entry.pack(side="left", padx=2)

        def clear_range():
            x_min_var.set("")
            x_max_var.set("")
            redraw()

        range_button_frame = ctk.CTkFrame(side_panel, fg_color="transparent")
        range_button_frame.pack(fill="x", padx=10, pady=(4, 0))
        ctk.CTkButton(range_button_frame, text="Apply Range", height=26,
                      command=redraw).pack(side="left", expand=True, fill="x", padx=2)
        ctk.CTkButton(range_button_frame, text="Full Range", height=26,
                      fg_color="#555555", hover_color="#666666",
                      command=clear_range).pack(side="left", expand=True, fill="x", padx=2)

        for entry in (min_entry, max_entry):
            entry.bind("<Return>", lambda _event: redraw())

        status_widget["label"] = ctk.CTkLabel(side_panel, textvariable=range_status_var,
                                              font=("Arial", 10), text_color="#aaaaaa")
        status_widget["label"].pack(padx=10, pady=(2, 0))
        ctk.CTkLabel(side_panel, text="Drag to pan • Scroll to zoom • Double-click to reset",
                     font=("Arial", 10), text_color="#aaaaaa").pack(padx=10, pady=(0, 4))

        ctk.CTkLabel(side_panel, text="Filters", font=("Arial", 12, "bold")).pack(pady=(8, 0))
        filter_frame = ctk.CTkFrame(side_panel, fg_color="transparent")
        filter_frame.pack(fill="x", padx=10, pady=(2, 0))
        ctk.CTkLabel(filter_frame, text="Bandwidth", font=("Arial", 11), width=64).pack(side="left")
        bandwidth_entry = ctk.CTkEntry(filter_frame, textvariable=bandwidth_var, width=70)
        bandwidth_entry.pack(side="left", padx=2)
        ctk.CTkButton(filter_frame, text="Add Filter", height=26,
                      command=lambda: add_filter()).pack(side="left", expand=True, fill="x", padx=2)

        def add_filter():
            """Highlights ±bandwidth/2 around the point last clicked on the graph."""
            if selection["wavenumber"] is None:
                set_filter_status("Click a point on the graph first", error=True)
                return

            text = bandwidth_var.get().strip().replace(",", ".")
            try:
                bandwidth = float(text)
            except ValueError:
                set_filter_status(f"Bandwidth: '{text}' is not a number", error=True)
                return
            if bandwidth <= 0:
                set_filter_status("Bandwidth must be greater than zero", error=True)
                return

            use_wavenumber = x_mode_var.get().startswith("Wavenumber")
            unit = "cm^-1" if use_wavenumber else "nm"
            centre_display = selection["display_x"]
            low_display = centre_display - bandwidth / 2.0
            high_display = centre_display + bandwidth / 2.0

            if use_wavenumber:
                low_cm, high_cm = low_display, high_display
            else:
                # In nm the band edges invert when converted to wavenumber.
                if low_display <= 0:
                    set_filter_status("Bandwidth reaches past 0 nm", error=True)
                    return
                low_cm, high_cm = 1.0e7 / high_display, 1.0e7 / low_display

            filter_bands.append({
                "low_cm": low_cm,
                "high_cm": high_cm,
                "centre_cm": selection["wavenumber"],
                "label": f"{centre_display:.1f} ± {bandwidth / 2:.1f} {unit}",
            })
            set_filter_status(f"{len(filter_bands)} filter(s): "
                              f"{low_display:.1f} - {high_display:.1f} {unit}")
            redraw()

        def clear_filters():
            filter_bands.clear()
            set_filter_status("Click a point on the graph, then Add Filter")
            redraw()

        ctk.CTkButton(side_panel, text="Clear Filters", height=24, font=("Arial", 11),
                      fg_color="#555555", hover_color="#666666",
                      command=clear_filters).pack(fill="x", padx=10, pady=(4, 0))
        bandwidth_entry.bind("<Return>", lambda _event: add_filter())

        filter_status_widget["label"] = ctk.CTkLabel(side_panel, textvariable=filter_status_var,
                                                     font=("Arial", 10), text_color="#aaaaaa",
                                                     wraplength=340)
        filter_status_widget["label"].pack(padx=10, pady=(2, 0))

        def set_all(visible):
            for sample in samples:
                sample["visible_var"].set(visible)
            redraw()

        bulk_frame = ctk.CTkFrame(side_panel, fg_color="transparent")
        bulk_frame.pack(fill="x", padx=10, pady=(5, 0))
        ctk.CTkButton(bulk_frame, text="Show All", width=80,
                      command=lambda: set_all(True)).pack(side="left", expand=True, padx=2)
        ctk.CTkButton(bulk_frame, text="Hide All", width=80, fg_color="#555555", hover_color="#666666",
                      command=lambda: set_all(False)).pack(side="left", expand=True, padx=2)
        ctk.CTkButton(bulk_frame, text="Reset Zoom", width=90, fg_color="#555555", hover_color="#666666",
                      command=reset_zoom).pack(side="left", expand=True, padx=2)

        ctk.CTkButton(side_panel, text="Load Text File(s)...", height=28,
                      fg_color="#a15c00", hover_color="#c47200",
                      command=lambda: load_external_files()).pack(fill="x", padx=10, pady=(6, 0))

        list_frame = ctk.CTkScrollableFrame(side_panel)
        list_frame.pack(fill="both", expand=True, padx=10, pady=10)

        def add_sample_row(sample):
            row = ctk.CTkFrame(list_frame, fg_color="transparent")
            row.pack(fill="x", pady=2)
            sample["row"] = row
            ctk.CTkLabel(row, text="■", text_color=sample["color"],
                         font=("Arial", 16), width=18).pack(side="left")
            # The Info button claims its space before the label so a long
            # composition string can never push it out of the row.
            ctk.CTkButton(
                row, text="Info", width=48, height=24, font=("Arial", 11),
                fg_color="#3a3a3a", hover_color="#4a4a4a",
                command=lambda s=sample, h=reference_headers: self.show_measured_sample_details(s, h),
            ).pack(side="right", padx=(4, 0))
            short_label = sample["label"]
            if len(short_label) > 34:
                short_label = short_label[:33] + "…"
            ctk.CTkCheckBox(
                row,
                text=short_label,
                variable=sample["visible_var"],
                command=redraw,
                font=("Arial", 11),
                checkbox_width=18,
                checkbox_height=18,
            ).pack(side="left", padx=(2, 4), fill="x", expand=True)

        for sample in samples:
            add_sample_row(sample)

        def load_external_files():
            """Loads arbitrary logs in the same format and adds them to the list."""
            popup.attributes("-topmost", False)  # Otherwise the dialog opens behind.
            try:
                paths = filedialog.askopenfilenames(
                    parent=popup,
                    title="Load spectrum log(s)",
                    initialdir=os.path.abspath(LOG_DIR) if os.path.isdir(LOG_DIR) else os.getcwd(),
                    filetypes=[("Spectrum logs", "*.txt"), ("All files", "*.*")],
                )
            finally:
                popup.attributes("-topmost", True)
            if not paths:
                return

            fallback = "absorbance" if y_mode_var.get() == "Absorbance" else "transmittance"
            added = 0
            problems = []
            for path in paths:
                wavenumber, _values, header = self.get_measured_spectrum(path)
                if wavenumber.size == 0:
                    problems.append(f"{os.path.basename(path)}: no data rows")
                    continue

                mode = infer_mode(path, header, default=fallback)
                key = external_group_key(path)

                existing = next((s for s in samples if s.get("external_key") == key), None)
                if existing is not None:
                    # Second half of an abs/trans pair, or a reload.
                    existing["paths"][mode] = path
                    added += 1
                    continue

                index = len(samples)
                sample = {
                    "number": None,
                    "variant": "",
                    "paths": {"absorbance": None, "transmittance": None},
                    "reference": {},
                    "external_key": key,
                    "source_path": path,
                    "key": strip_mode_tokens(path),
                    "label": strip_mode_tokens(path),
                    "color": MEASURED_PALETTE[index % len(MEASURED_PALETTE)],
                    "linestyle": MEASURED_LINESTYLES[(index // len(MEASURED_PALETTE)) % len(MEASURED_LINESTYLES)],
                    "visible_var": ctk.BooleanVar(value=True),
                }
                sample["paths"][mode] = path
                samples.append(sample)
                add_sample_row(sample)
                added += 1

            if added:
                redraw()
            if problems:
                set_status("; ".join(problems), error=True)
            elif added:
                set_status(f"Loaded {added} file(s)")

        ctk.CTkLabel(
            side_panel,
            text=("Click a trace to read off a value.\n"
                  "Scroll to zoom (Ctrl = X only, Shift = Y only).\n"
                  "Double-click or Reset Zoom to go back.\n"
                  "X range clips the data; leave a box blank for open-ended.\n"
                  "Loaded files pair up by name (_abs with _trans).\n"
                  "Add Filter highlights ±bandwidth/2 around the clicked point.\n"
                  "Repeated colours differ by line style."),
            font=("Arial", 10), text_color="#aaaaaa", justify="left",
        ).pack(padx=10, pady=(0, 10), anchor="w")

        redraw()
        canvas.draw()

    def show_measured_sample_details(self, sample, reference_headers):
        """Shows the reference-sheet row and log details for one sample."""
        popup = ctk.CTkToplevel(self)
        popup.title(f"Details - {sample['label']}")
        popup.geometry("640x520")
        popup.attributes("-topmost", True)

        container = ctk.CTkScrollableFrame(popup)
        container.pack(fill="both", expand=True, padx=15, pady=15)
        container.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(container, text=sample["label"], font=("Arial", 16, "bold"),
                     text_color=sample["color"], wraplength=560, justify="left").grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 10))

        row_index = 1

        def add_row(key, value):
            nonlocal row_index
            ctk.CTkLabel(container, text=key, font=("Arial", 12, "bold"),
                         anchor="w").grid(row=row_index, column=0, sticky="nw", padx=(0, 10), pady=3)
            ctk.CTkLabel(container, text=value, font=("Arial", 12), anchor="w",
                         wraplength=380, justify="left").grid(row=row_index, column=1, sticky="w", pady=3)
            row_index += 1

        reference = sample.get("reference") or {}
        if reference:
            ctk.CTkLabel(container, text="Reference sheet (sample_references.xlsx)",
                         font=("Arial", 13, "bold"), text_color="#aaaaaa").grid(
                row=row_index, column=0, columnspan=2, sticky="w", pady=(10, 4))
            row_index += 1
            for _column, header in reference_headers:
                add_row(header, reference.get(header, "-") or "-")
        else:
            add_row("Reference sheet", "No matching row found.")

        ctk.CTkLabel(container, text="Measured logs", font=("Arial", 13, "bold"),
                     text_color="#aaaaaa").grid(row=row_index, column=0, columnspan=2,
                                                sticky="w", pady=(14, 4))
        row_index += 1

        if sample.get("source_path"):
            add_row("Loaded from", os.path.dirname(os.path.abspath(sample["source_path"])))
        if sample["variant"]:
            add_row("Measurement variant", sample["variant"].replace("_", " "))
        add_row("Line style", MEASURED_LINESTYLE_GLYPHS.get(sample["linestyle"], sample["linestyle"]))

        for mode_key, mode_label in (("absorbance", "Absorbance log"), ("transmittance", "Transmittance log")):
            path = sample["paths"].get(mode_key)
            if not path:
                add_row(mode_label, "not available")
                continue
            wavenumber, values, header = self.get_measured_spectrum(path)
            add_row(mode_label, os.path.basename(path))
            if wavenumber.size:
                add_row(
                    "   range",
                    f"{wavenumber.min():.1f} - {wavenumber.max():.1f} cm^-1"
                    f"  ({wavenumber_to_nm(wavenumber.max()):.0f} - {wavenumber_to_nm(wavenumber.min()):.0f} nm)",
                )
                add_row("   points", f"{wavenumber.size}")
                add_row("   y range", f"{values.min():.4g} to {values.max():.4g} "
                                      f"({header.get('YUNITS', 'unknown units')})")

    def show_thickness_solver(self):
        channels = self.get_channel_definitions()
        materials = self.get_unique_stack_materials()

        if not channels or not materials:
            popup = ctk.CTkToplevel(self)
            popup.title("Solve Thicknesses")
            popup.geometry("560x180")
            popup.attributes("-topmost", True)
            ctk.CTkLabel(
                popup,
                text="Add at least one channel and one non-air material layer before solving.",
                font=("Arial", 14, "bold"),
            ).pack(expand=True, padx=20, pady=20)
            return

        matrix, reference_signals = self.build_effective_channel_matrix(channels, materials)
        source_obj = LightSource("Source", self.wl, self.get_source_spectra())
        web_stack_objs = self.build_web_stack_objects()

        measured_signals = []
        for channel in channels:
            filter_spectra = gaussian_bandpass(self.wl, channel["center"], channel["width"], 1.0)
            sensor_spectra = self.get_sensor_spectra(channel["sensor_type"])
            filter_obj = OpticalFilter(f"Filter_{channel['center']}", self.wl, filter_spectra)
            sensor_obj = Sensor("Sensor", self.wl, sensor_spectra)
            results = run_simulation(self.wl, source_obj, web_stack_objs, filter_obj, sensor_obj)
            measured_signals.append(results["final_signal"])

        measured_signals = np.asarray(measured_signals, dtype=float)
        safe_reference = np.clip(reference_signals, 1e-12, None)
        safe_measured = np.clip(measured_signals, 1e-12, None)
        absorbance = -np.log(safe_measured / safe_reference)

        active_rows = reference_signals > 1e-9
        solve_matrix = matrix[active_rows, :]
        solve_absorbance = absorbance[active_rows]

        rank = int(np.linalg.matrix_rank(solve_matrix, tol=1e-9)) if solve_matrix.size else 0
        singular_values = np.linalg.svd(solve_matrix, compute_uv=False) if solve_matrix.size else np.array([])
        nonzero_singular_values = singular_values[singular_values > 1e-12]
        condition = (
            nonzero_singular_values[0] / nonzero_singular_values[-1]
            if len(nonzero_singular_values) >= 2
            else float("inf")
        )

        if solve_matrix.size:
            least_squares_solution, residuals, _, _ = np.linalg.lstsq(solve_matrix, solve_absorbance, rcond=None)
            nnls_solution, nnls_residual_norm = nnls(solve_matrix, solve_absorbance)
            predicted_absorbance = solve_matrix @ nnls_solution
            residual_vector = solve_absorbance - predicted_absorbance
            rms_residual = float(np.sqrt(np.mean(residual_vector ** 2))) if residual_vector.size else 0.0
        else:
            least_squares_solution = np.zeros(len(materials))
            nnls_solution = np.zeros(len(materials))
            nnls_residual_norm = float("nan")
            residuals = np.array([])
            rms_residual = float("nan")

        actual_totals = {material: 0.0 for material in materials}
        for layer in self.web_layers:
            material = layer["mat_var"].get()
            if material not in actual_totals:
                continue
            try:
                actual_totals[material] += float(layer["thick_var"].get())
            except ValueError:
                pass

        popup = ctk.CTkToplevel(self)
        popup.title("Solve Thicknesses")
        popup.geometry("1000x720")
        popup.attributes("-topmost", True)

        ctk.CTkLabel(
            popup,
            text="Solved material thicknesses from selected channels and simulated signals",
            font=("Arial", 16, "bold"),
        ).pack(pady=(10, 4))

        condition_text = f"{condition:.3g}" if np.isfinite(condition) else "inf"
        summary = (
            f"Channels used: {int(np.count_nonzero(active_rows))}/{len(channels)}    "
            f"Materials: {len(materials)}    Rank: {rank}/{len(materials)}    "
            f"Condition: {condition_text}"
        )
        ctk.CTkLabel(popup, text=summary).pack(pady=(0, 8))

        text_box = ctk.CTkTextbox(popup, font=("Consolas", 12))
        text_box.pack(fill="both", expand=True, padx=10, pady=(5, 10))

        lines = [
            "This solves A = K * thickness, where A = -ln(measured signal / open-beam signal).",
            "NNLS is the recommended estimate because thickness cannot be negative.",
            "Interface/Fresnel effects are included in the simulated measured signals, so large residuals can indicate model mismatch or weak channels.",
            "",
            "Material".ljust(24) + "Actual mm".rjust(14) + "NNLS mm".rjust(14) + "LeastSq mm".rjust(14),
        ]

        for idx, material in enumerate(materials):
            lines.append(
                material[:23].ljust(24)
                + f"{actual_totals[material]:14.6g}"
                + f"{nnls_solution[idx]:14.6g}"
                + f"{least_squares_solution[idx]:14.6g}"
            )

        ls_residual_text = f"{float(np.sum(residuals)):.6g}" if residuals.size else "n/a"
        lines.extend([
            "",
            f"NNLS residual norm: {nnls_residual_norm:.6g}",
            f"NNLS RMS absorbance residual: {rms_residual:.6g}",
            f"Least-squares summed residual: {ls_residual_text}",
            "",
            "Channel".ljust(30) + "Signal".rjust(14) + "Open beam".rjust(14) + "Absorbance".rjust(14),
        ])

        for idx, channel in enumerate(channels):
            active_marker = "" if active_rows[idx] else "  low-ref"
            name = f"{channel['name']} {channel['center']:.0f}nm"
            lines.append(
                name[:29].ljust(30)
                + f"{measured_signals[idx]:14.6g}"
                + f"{reference_signals[idx]:14.6g}"
                + f"{absorbance[idx]:14.6g}"
                + active_marker
            )

        if len(channels) < len(materials):
            lines.extend(["", "Warning: fewer channels than materials. The solve is underdetermined."])
        elif rank < len(materials):
            lines.extend(["", "Warning: active channel matrix is rank deficient. Some materials cannot be separated."])
        elif condition > 100:
            lines.extend(["", "Warning: high condition number. Small signal errors can create large thickness errors."])

        text_box.insert("1.0", "\n".join(lines))
        text_box.configure(state="disabled")

    def show_channel_matrix(self):
        channels = self.get_channel_definitions()
        materials = self.get_unique_stack_materials()

        if not channels or not materials:
            popup = ctk.CTkToplevel(self)
            popup.title("Channel Matrix")
            popup.geometry("560x180")
            popup.attributes("-topmost", True)
            ctk.CTkLabel(
                popup,
                text="Add at least one channel and one non-air truth-stack material before viewing the matrix.",
                font=("Arial", 14, "bold"),
            ).pack(expand=True, padx=20, pady=20)
            return

        matrix, channel_weights = self.build_effective_channel_matrix(channels, materials)
        singular_values = np.linalg.svd(matrix, compute_uv=False)
        nonzero_singular_values = singular_values[singular_values > 1e-12]
        rank = int(np.linalg.matrix_rank(matrix, tol=1e-9))
        condition = nonzero_singular_values[0] / nonzero_singular_values[-1] if len(nonzero_singular_values) >= 2 else float("inf")

        popup = ctk.CTkToplevel(self)
        popup.title("Channel Matrix")
        popup.geometry("1100x850")
        popup.attributes("-topmost", True)

        header = ctk.CTkLabel(popup, text="Effective alpha matrix for truth-stack materials", font=("Arial", 16, "bold"))
        header.pack(pady=(10, 4))

        summary_text = f"Channels: {len(channels)}    Truth-stack materials: {len(materials)}    Rank: {rank}    Condition: {condition:.3g}"
        ctk.CTkLabel(popup, text=summary_text).pack(pady=(0, 8))
        ctk.CTkLabel(
            popup,
            text="Curve-display checkboxes only affect spectra plots; this matrix only includes materials currently used in the web stack.",
        ).pack(pady=(0, 8))

        # --- Filter Selection UI ---
        combo_container = ctk.CTkFrame(popup, fg_color="#333333")
        combo_container.pack(fill="x", padx=10, pady=5)

        input_row = ctk.CTkFrame(combo_container, fg_color="transparent")
        input_row.pack(fill="x", padx=5, pady=5)

        ctk.CTkLabel(input_row, text="Evaluate Best Combinations for N Filters:").pack(side="left", padx=5)
        n_var = ctk.StringVar(value=str(min(len(materials), len(channels))))
        ctk.CTkEntry(input_row, textvariable=n_var, width=50).pack(side="left", padx=5)

        self.top_combos = []
        combo_choice_var = ctk.IntVar(value=0)

        results_row = ctk.CTkFrame(combo_container, fg_color="transparent")
        results_row.pack(fill="x", padx=5, pady=5)

        def calculate_best_combo():
            for widget in results_row.winfo_children(): widget.destroy()

            try: n = int(n_var.get())
            except ValueError: return
            if n <= 0 or n > len(channels): return

            combo_results = []
            for combo in itertools.combinations(range(len(channels)), n):
                sub_mat = matrix[list(combo), :]
                sub_rank = int(np.linalg.matrix_rank(sub_mat, tol=1e-9))

                if sub_rank < len(materials):
                    cond = float("inf")
                else:
                    s_vals = np.linalg.svd(sub_mat, compute_uv=False)
                    nonzero_s = s_vals[s_vals > 1e-12]
                    cond = nonzero_s[0] / nonzero_s[-1] if len(nonzero_s) >= 2 else float("inf")

                combo_results.append((cond, combo))

            combo_results.sort(key=lambda x: x[0])
            self.top_combos = combo_results[:5]
            combo_choice_var.set(0)

            ctk.CTkLabel(results_row, text="Top 5 Filter Combinations:", font=("Arial", 12, "bold")).pack(anchor="w", padx=5, pady=(0, 5))

            for i, (cond, combo) in enumerate(self.top_combos):
                ch_names = [f"{channels[idx]['center']:.0f}nm" for idx in combo]
                total_price = sum([channels[idx]['price'] for idx in combo])

                lbl_text = f"Rank {i+1} (Condition: {cond:.2f} | Price: ${total_price:.2f}):   {', '.join(ch_names)}"
                rb = ctk.CTkRadioButton(results_row, text=lbl_text, variable=combo_choice_var, value=i)
                rb.pack(anchor="w", padx=10, pady=2)

        def apply_combo():
            if not self.top_combos: return
            choice_idx = combo_choice_var.get()
            best_combo_indices = self.top_combos[choice_idx][1]
            selected_channels = [channels[idx] for idx in best_combo_indices]

            for ch in self.sensor_channels: ch["frame"].destroy()
            self.sensor_channels.clear()

            for ch in selected_channels:
                self.add_sensor_ui(center=str(ch["center"]), width=str(ch["width"]), name=ch["name"], sensor_type=ch["sensor_type"], price=str(ch["price"]))

            self.run_live_simulation()
            popup.destroy()

        ctk.CTkButton(input_row, text="Calculate Options", command=calculate_best_combo).pack(side="left", padx=10)
        ctk.CTkButton(input_row, text="Apply Selected Combination", command=apply_combo, fg_color="#00aa00", hover_color="#008800").pack(side="left", padx=10)

        fig_width = max(8.0, 1.25 * len(materials) + 3.0)
        fig_height = max(4.8, 0.52 * len(channels) + 2.6)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor='#2b2b2b')
        ax.set_facecolor('#2b2b2b')

        if np.nanmax(matrix) > 0:
            display_matrix = np.log10(matrix + 1e-9)
            image_label = "log10(effective alpha + 1e-9)"
        else:
            display_matrix = matrix
            image_label = "effective alpha"

        im = ax.imshow(display_matrix, aspect="auto", cmap="viridis")
        ax.set_title(image_label, color="white")
        ax.set_xticks(np.arange(len(materials)))
        ax.set_yticks(np.arange(len(channels)))
        ax.set_xticklabels(materials, color="white", rotation=35 if len(materials) > 5 else 0, ha="right" if len(materials) > 5 else "center")
        ax.set_yticklabels([f"{ch['name']} ({ch['center']:.0f} nm)" for ch in channels], color="white")
        ax.tick_params(colors="white")

        for row_idx in range(len(channels)):
            for col_idx in range(len(materials)):
                value = matrix[row_idx, col_idx]
                ax.text(col_idx, row_idx, f"{value:.3g}", ha="center", va="center", color="white", fontsize=8)

        fig.colorbar(im, ax=ax).ax.tick_params(colors="white")
        fig.tight_layout()

        plot_container = ctk.CTkFrame(popup)
        plot_container.pack(fill="both", expand=True, padx=10, pady=5)

        scroll_canvas = tk.Canvas(plot_container, bg="#2b2b2b", highlightthickness=0)
        y_scroll = tk.Scrollbar(plot_container, orient="vertical", command=scroll_canvas.yview)
        x_scroll = tk.Scrollbar(plot_container, orient="horizontal", command=scroll_canvas.xview)
        scroll_canvas.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)

        y_scroll.pack(side="right", fill="y")
        x_scroll.pack(side="bottom", fill="x")
        scroll_canvas.pack(side="left", fill="both", expand=True)

        plot_frame = tk.Frame(scroll_canvas, bg="#2b2b2b")
        scroll_canvas.create_window((0, 0), window=plot_frame, anchor="nw")

        canvas = FigureCanvasTkAgg(fig, master=plot_frame)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack()
        canvas_widget.configure(width=int(fig_width * fig.dpi), height=int(fig_height * fig.dpi))

        def update_scroll_region(_event=None):
            scroll_canvas.configure(scrollregion=scroll_canvas.bbox("all"))

        plot_frame.bind("<Configure>", update_scroll_region)
        canvas.draw()
        update_scroll_region()

        table_lines = [
            "Effective alpha values are in mm^-1.",
            "Use absorbance A = -ln(signal/reference). Then approximately A = K * thickness.",
            "",
            "Channel".ljust(28) + "".join(material.rjust(14) for material in materials) + "    Weight"
        ]
        for row_idx, channel in enumerate(channels):
            name = f"{channel['name']} {channel['center']:.0f}nm"
            values = "".join(f"{matrix[row_idx, col_idx]:14.5g}" for col_idx in range(len(materials)))
            table_lines.append(name[:27].ljust(28) + values + f"    {channel_weights[row_idx]:.3g}")

        table_lines.append("")
        if len(channels) < len(materials): table_lines.append("Warning: fewer channels than materials.")
        elif rank < len(materials): table_lines.append("Warning: matrix is rank deficient.")
        elif condition > 100: table_lines.append("Warning: high condition number.")
        else: table_lines.append("Matrix looks reasonably separable for these modeled curves.")

        text_box = ctk.CTkTextbox(popup, height=170, font=("Consolas", 12))
        text_box.pack(fill="x", padx=10, pady=(5, 10))
        text_box.insert("1.0", "\n".join(table_lines))
        text_box.configure(state="disabled")

    def show_ranked_combinations(self):
        """Ranks source/filter/sensor combinations from the component database."""
        def price_value(component, component_kind):
            value = component.get("price_usd")
            if value is None:
                defaults = self.component_database.get("price_defaults_usd", {})
                value = defaults.get(component_kind, {}).get(component.get("manufacturer"))
                if value is None:
                    return None
            try:
                return float(value)
            except (TypeError, ValueError):
                return None

        def price_label(component, component_kind):
            value = price_value(component, component_kind)
            if value is None:
                return "price: quote/unknown"
            label = f"price: ${value:,.2f}"
            status = component.get("price_status", "default estimate")
            if status:
                label += f" ({status})"
            return label

        def combo_price_summary(result):
            components = [("source", result["source"])]
            for channel in result["channels"]:
                components.append(("filter", channel["filter"]))
                components.append(("sensor", channel["sensor"]))

            known_total = 0.0
            unknown_count = 0
            for component_kind, component in components:
                value = price_value(component, component_kind)
                if value is None:
                    unknown_count += 1
                else:
                    known_total += value

            if unknown_count:
                return f"${known_total:,.2f} + {unknown_count} quote/unknown item(s)"
            return f"${known_total:,.2f}"

        materials = [
            material for material in material_names_from_stack(self.material_library, self.web_layers)
            if self.is_material_displayed(material)
        ]
        if not materials:
            popup = ctk.CTkToplevel(self)
            popup.title("Ranked Filter / Source / Sensor Combos")
            popup.geometry("520x160")
            popup.attributes("-topmost", True)
            ctk.CTkLabel(
                popup,
                text="Select at least one displayed material before ranking combinations.",
                font=("Arial", 14, "bold"),
            ).pack(expand=True, padx=20, pady=20)
            return

        channel_count = max(len(materials), 1)
        results = rank_orthogonal_combinations(
            self.wl,
            self.material_library,
            self.component_database,
            materials=materials,
            channel_count=channel_count,
            top_n=25,
            beam_width=300,
        )

        popup = ctk.CTkToplevel(self)
        popup.title("Ranked Filter / Source / Sensor Combos")
        popup.geometry("1150x760")
        popup.attributes("-topmost", True)

        ctk.CTkLabel(
            popup,
            text="Ranked database combinations for the current material stack",
            font=("Arial", 16, "bold"),
        ).pack(pady=(10, 4))

        summary = (
            f"Materials: {', '.join(materials)}    "
            f"Channels per combo: {channel_count}    "
            f"Sources: {len(self.component_database.get('sources', []))}    "
            f"Filters: {len(self.component_database.get('filters', []))}    "
            f"Sensors: {len(self.component_database.get('sensors', []))}"
        )
        ctk.CTkLabel(popup, text=summary).pack(pady=(0, 8))

        text_box = ctk.CTkTextbox(popup, font=("Consolas", 12))
        text_box.pack(fill="both", expand=True, padx=10, pady=(5, 10))

        if not results:
            text_box.insert("1.0", "No viable combinations found. Check source/filter/sensor wavelength overlap.")
            text_box.configure(state="disabled")
            return

        lines = []
        lines.append("Score favors full-rank matrices, low material-column similarity, good conditioning, and live detector signal.")
        lines.append("Use these as modeled shortlists; final hardware choices still need calibration and signal/noise checks.")
        lines.append("")

        for idx, result in enumerate(results, start=1):
            condition = result["condition"]
            condition_text = f"{condition:.3g}" if np.isfinite(condition) else "inf"
            lines.append(
                f"{idx:02d}. Score {result['score']:.1f} | Rank {result['rank']}/{len(materials)} | "
                f"Condition {condition_text} | Orthogonality {result['orthogonality']:.3f} | "
                f"Set price {combo_price_summary(result)}"
            )
            source_url = result["source"].get("url", "")
            lines.append(f"    Source: {result['source']['name']} | {price_label(result['source'], 'source')}")
            if source_url:
                lines.append(f"       link: {source_url}")
            for channel_idx, channel in enumerate(result["channels"], start=1):
                filter_def = channel["filter"]
                sensor_def = channel["sensor"]
                lines.append(
                    f"    Ch {channel_idx}: {filter_def['name']} "
                    f"({float(filter_def['center_nm']):.0f} nm / {float(filter_def['fwhm_nm']):.0f} nm) "
                    f"+ {sensor_def['name']} | weight {result['weights'][channel_idx - 1]:.3g}"
                )
                lines.append(f"       filter {price_label(filter_def, 'filter')} | link: {filter_def.get('url', 'n/a')}")
                lines.append(f"       sensor {price_label(sensor_def, 'sensor')} | link: {sensor_def.get('url', 'n/a')}")
            lines.append("    Effective alpha matrix rows are channels, columns are " + ", ".join(materials))
            for row in result["matrix"]:
                lines.append("       " + " ".join(f"{value:10.4g}" for value in row))
            lines.append("")

        text_box.insert("1.0", "\n".join(lines))
        text_box.configure(state="disabled")

    # --- Core Physics Engine ---
    def run_live_simulation(self):
        source_spectra = self.get_source_spectra()
        source_obj = LightSource("Source", self.wl, source_spectra)
        web_stack_objs = self.build_web_stack_objects()

        for ax in (self.ax_top, self.ax_bot):
            ax.clear()
            ax.set_facecolor('#2b2b2b')
            ax.tick_params(colors='white')
            for spine in ax.spines.values(): spine.set_color('gray')

        self.ax_top.set_title("System Optics: Source Emittance & Web Transmission", color='white')
        self.ax_top.set_ylabel("Transmission / Intensity", color='white')
        self.ax_bot.set_title("Filtered Signals Reaching Detectors", color='white')
        self.ax_bot.set_xlabel("Wavelength (nm)", color='white')

        self.ax_top.plot(self.wl, source_spectra, color='#aaaaaa', linestyle="--", label="Source Output")

        first_channel = True

        for ch in self.sensor_channels:
            try: c_wl, fwhm = float(ch["center_var"].get()), float(ch["width_var"].get())
            except ValueError: c_wl, fwhm = 2000, 10

            filter_spectra = gaussian_bandpass(self.wl, c_wl, fwhm, 1.0)
            sensor_spectra = self.get_sensor_spectra(ch["sensor_type_var"].get())

            filter_obj = OpticalFilter(f"Filter_{c_wl}", self.wl, filter_spectra)
            sensor_obj = Sensor("Sensor", self.wl, sensor_spectra)

            results = run_simulation(self.wl, source_obj, web_stack_objs, filter_obj, sensor_obj)
            final_signal, channel_spectra = results["final_signal"], results["spectra"]["signal_spectrum"]

            if first_channel:
                T_total = results["spectra"]["bulk_transmission"] * results["spectra"]["interface_transmission"]
                self.ax_top.plot(self.wl, T_total, color='white', label="Transmitted Spectrum (Bulk + Fresnel)")
                first_channel = False

            ch["lbl_readout"].configure(text=f"Signal: {final_signal:.2f}")

            name, color = ch["name_var"].get(), ch["color"]
            self.ax_top.fill_between(self.wl, 0, filter_spectra, color=color, alpha=0.2)
            self.ax_bot.plot(self.wl, channel_spectra, color=color, label=f"{name} (Signal: {final_signal:.1f})")
            self.ax_bot.fill_between(self.wl, 0, channel_spectra, color=color, alpha=0.5)

        self.ax_top.legend(facecolor='#333333', edgecolor='white', labelcolor='white', loc='upper right')
        if self.sensor_channels:
            self.ax_bot.legend(facecolor='#333333', edgecolor='white', labelcolor='white', loc='upper right')

        self.fig.tight_layout()
        self.canvas.draw()

if __name__ == "__main__":
    app = WebGaugingApp()
    app.mainloop()
