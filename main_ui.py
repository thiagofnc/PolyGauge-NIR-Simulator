import customtkinter as ctk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import os 

# Import the backend components
from Components import LightSource, OpticalFilter, Sensor, MaterialLayer
from Simulation import run_simulation

# --- Physics & Data Helpers ---
def gaussian_peak(wl, center, width, amplitude):
    return amplitude * np.exp(-((wl - center) ** 2) / (2 * width ** 2))

def gaussian_bandpass(wl, center, fwhm, peak=1.0):
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return gaussian_peak(wl, center, sigma, peak)

def band_model(wl, bands):
    """Build a simple absorption model from NIR band centers.

    Amplitudes are engineering-scale alpha values in mm^-1. They are intended
    for comparing band placement and channel sensitivity, not for calibration.
    """
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
    """
    Parses standard YAML-like optical database files (e.g., refractiveindex.info).
    Extracts Wavelength, n, and k, and interpolates them to the master array.
    """
    if not os.path.exists(filepath):
        print(f"Warning: Could not find '{filepath}'.")
        return None, None
        
    raw_wl, raw_n, raw_k = [], [], []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split()
            # We are looking for lines with at least 2 numbers (wl and n)
            if len(parts) >= 2:
                try:
                    # If this succeeds, we are looking at real data, not a header
                    wl = float(parts[0])
                    n = float(parts[1])
                    raw_wl.append(wl)
                    raw_n.append(n)
                    
                    # If there is a 3rd column, it is 'k'
                    if len(parts) >= 3:
                        raw_k.append(float(parts[2]))
                except ValueError:
                    # It hit text like "REFERENCES:" or "data:". Skip it.
                    pass
                    
    if not raw_wl:
        return None, None

    # Convert wavelengths from micrometers to nanometers
    raw_wl = np.array(raw_wl) * 1000.0
    raw_n = np.array(raw_n)
    
    # Interpolate n onto master wavelengths
    interp_n = np.interp(master_wl, raw_wl, raw_n, left=raw_n[0], right=raw_n[-1])
    
    # Interpolate k if it exists
    if raw_k:
        raw_k = np.array(raw_k)
        interp_k = np.interp(master_wl, raw_wl, raw_k, left=raw_k[0], right=raw_k[-1])
    else:
        interp_k = None
        
    return interp_n, interp_k

# Set UI Theme
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

class WebGaugingApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("PolyGauge NIR Simulator")
        self.geometry("1400x850")

        # Master Data
        self.wl = np.arange(1000, 4001, 1.0)
        self.web_layers = []
        self.sensor_channels = []
        
        self.colors = ['#00ffcc', '#ff3366', '#ffcc00', '#cc33ff', '#33ccff']


        #   --- Pre-loaded Material Library ---
        self.material_library = {
            "Air": {"alpha": np.zeros_like(self.wl), "n": 1.0},
            # NIR engineering band models for absorbance-band comparison.
            # See README "Material Data Notes" for source references.
            "PE": {
                "alpha": band_model(self.wl, [
                    (1210, 35, 0.20),
                    (1730, 24, 0.85),
                    (1764, 24, 0.55),
                    (2310, 32, 2.00),
                    (2350, 32, 1.35),
                ]),
                "n": 1.51,
            },
            "EVOH": {
                "alpha": band_model(self.wl, [
                    (1410, 28, 0.45),
                    (2012, 38, 0.90),
                    (2092, 48, 1.80),
                    (2310, 35, 0.45),
                    (3000, 90, 18.0),
                    (3305, 110, 12.0),
                ]),
                "n": 1.52,
            },
            "Nylon": {
                "alpha": band_model(self.wl, [
                    (1535, 30, 1.20),
                    (2040, 42, 1.50),
                    (2300, 35, 0.55),
                    (2355, 35, 0.55),
                    (3030, 75, 18.0),
                    (3410, 80, 8.0),
                    (3500, 85, 6.0),
                ]),
                "n": 1.53,
            },
        }

        pe_n, pe_k = load_database_file("PE_data.yml", self.wl)
        if pe_k is not None:
            measured_region = self.wl >= 2600
            pe_alpha = self.material_library["PE"]["alpha"].copy()
            pe_alpha[measured_region] = alpha_from_k(self.wl, pe_k)[measured_region]

            pe_n_combined = np.full_like(self.wl, 1.51, dtype=float)
            if pe_n is not None:
                pe_n_combined[measured_region] = pe_n[measured_region]

            self.material_library["PE"] = {
                "alpha": pe_alpha,
                "n": pe_n_combined,
            }

        water_n, water_k = load_database_file("Water_data.yml", self.wl)
        if water_k is not None:
            self.material_library["Water"] = {
                "alpha": alpha_from_k(self.wl, water_k),
                "n": water_n if water_n is not None else 1.33,
            }
        else:
            self.material_library["Water"] = {
                "alpha": band_model(self.wl, [
                    (1450, 55, 12.0),
                    (1940, 70, 40.0),
                ]),
                "n": 1.33,
            }

        self.setup_ui()
        self.run_live_simulation()

    def setup_ui(self):
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=3)
        self.grid_rowconfigure(0, weight=1)

        # ==========================================
        # PANEL 1: SOURCE & SENSOR STACK (Left)
        # ==========================================
        self.left_panel = ctk.CTkFrame(self)
        self.left_panel.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        
        ctk.CTkLabel(self.left_panel, text="Light Source", font=("Arial", 18, "bold")).pack(pady=(10, 0))
        self.source_frame = ctk.CTkFrame(self.left_panel, fg_color="#2b2b2b")
        self.source_frame.pack(fill="x", padx=10, pady=5)
        
        self.source_type_var = ctk.StringVar(value="Blackbody (Halogen)")
        ctk.CTkOptionMenu(self.source_frame, variable=self.source_type_var, 
                          values=[
                              "Blackbody (Halogen)",
                              "MTE6114W-WRC LED (1460nm)",
                              "HPIR104 Thermal Emitter",
                              "Flat Emission (Ideal)",
                          ]).pack(pady=5)
        
        temp_frame = ctk.CTkFrame(self.source_frame, fg_color="transparent")
        temp_frame.pack(fill="x", pady=5)
        ctk.CTkLabel(temp_frame, text="Temp (K):").pack(side="left", padx=5)
        self.source_temp_var = ctk.StringVar(value="3000")
        ctk.CTkEntry(temp_frame, textvariable=self.source_temp_var, width=60).pack(side="left", padx=5)

        ctk.CTkLabel(self.left_panel, text="Sensor Channels", font=("Arial", 18, "bold")).pack(pady=(20, 0))
        self.sensors_container = ctk.CTkScrollableFrame(self.left_panel)
        self.sensors_container.pack(fill="both", expand=True, padx=5, pady=5)
        
        ctk.CTkButton(self.left_panel, text="+ Add Sensor Channel", command=self.add_sensor_ui).pack(pady=10)
        
        self.add_sensor_ui(center="2310", width="15", name="Measure (PE)")
        self.add_sensor_ui(center="2200", width="15", name="Reference")

        # ==========================================
        # PANEL 2: WEB STACK (Center)
        # ==========================================
        self.stack_panel = ctk.CTkFrame(self)
        self.stack_panel.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        
        ctk.CTkLabel(self.stack_panel, text="Web Material Stack", font=("Arial", 18, "bold")).pack(pady=10)
        self.layers_container = ctk.CTkScrollableFrame(self.stack_panel)
        self.layers_container.pack(fill="both", expand=True, padx=5, pady=5)
        
        ctk.CTkButton(self.stack_panel, text="+ Add Layer", command=self.add_layer_ui).pack(pady=10)
        ctk.CTkButton(self.stack_panel, text="Compare Absorbance Curves",
                      command=self.show_all_material_spectra).pack(pady=(0, 10))
        
        self.add_layer_ui(mat="PE", thick="0.05")
        self.add_layer_ui(mat="EVOH", thick="0.015")
        self.add_layer_ui(mat="PE", thick="0.05")

        # ==========================================
        # PANEL 3: GRAPH & GLOBAL CONTROLS (Right)
        # ==========================================
        self.data_panel = ctk.CTkFrame(self)
        self.data_panel.grid(row=0, column=2, padx=10, pady=10, sticky="nsew")
        
        self.fig, (self.ax_top, self.ax_bot) = plt.subplots(2, 1, figsize=(8, 6), facecolor='#2b2b2b', gridspec_kw={'height_ratios': [2, 1]})
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.data_panel)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

        ctk.CTkButton(self.data_panel, text="UPDATE SIMULATION", command=self.run_live_simulation, 
                      fg_color="#00aa00", hover_color="#008800", font=("Arial", 16, "bold"), height=40).pack(fill="x", padx=20, pady=10)
        ctk.CTkButton(self.data_panel, text="CHANNEL MATRIX", command=self.show_channel_matrix,
                      fg_color="#0055aa", hover_color="#0077cc", font=("Arial", 15, "bold"), height=36).pack(fill="x", padx=20, pady=(0, 10))

    # --- UI Builders ---

    def add_layer_ui(self, mat="PE", thick="0.05"):
        frame = ctk.CTkFrame(self.layers_container, fg_color="#333333")
        frame.pack(fill="x", pady=5, padx=5)
        
        mat_var = ctk.StringVar(value=mat)
        ctk.CTkOptionMenu(frame, variable=mat_var, values=list(self.material_library.keys()), width=90).pack(side="left", padx=5, pady=10)
        
        thick_var = ctk.StringVar(value=thick)
        ctk.CTkEntry(frame, textvariable=thick_var, width=60).pack(side="left", padx=5)
        ctk.CTkLabel(frame, text="mm").pack(side="left")
        
        # Delete Button (Red)
        ctk.CTkButton(frame, text="X", width=30, fg_color="#aa0000", hover_color="#ff0000",
                      command=lambda f=frame: self.remove_element(f, self.web_layers)).pack(side="right", padx=(5, 5))
        
        # View Spectral Data Button (Blue)
        ctk.CTkButton(frame, text="👁️", width=30, fg_color="#0055aa", hover_color="#0077ff",
                      command=lambda m=mat_var: self.show_material_spectra(m.get())).pack(side="right", padx=(5, 0))

        self.web_layers.append({"frame": frame, "mat_var": mat_var, "thick_var": thick_var})

    def add_sensor_ui(self, center="2310", width="15", name="Sensor"):
        color = self.colors[len(self.sensor_channels) % len(self.colors)]
        frame = ctk.CTkFrame(self.sensors_container, border_width=2, border_color=color, fg_color="#333333")
        frame.pack(fill="x", pady=5, padx=5)
        
        row1 = ctk.CTkFrame(frame, fg_color="transparent")
        row1.pack(fill="x", padx=5, pady=(5,0))
        name_var = ctk.StringVar(value=name)
        ctk.CTkEntry(row1, textvariable=name_var, width=100).pack(side="left")
        sensor_type_var = ctk.StringVar(value="InGaAs")
        ctk.CTkOptionMenu(row1, variable=sensor_type_var,
                          values=["InGaAs", "InAsSb (2.7-5.3um)", "MCT (MIR)", "Ideal (Flat)"],
                          width=130).pack(side="right")

        row2 = ctk.CTkFrame(frame, fg_color="transparent")
        row2.pack(fill="x", padx=5, pady=5)
        ctk.CTkLabel(row2, text="CWL:").pack(side="left")
        center_var = ctk.StringVar(value=center)
        ctk.CTkEntry(row2, textvariable=center_var, width=50).pack(side="left", padx=(0,10))
        ctk.CTkLabel(row2, text="FWHM:").pack(side="left")
        width_var = ctk.StringVar(value=width)
        ctk.CTkEntry(row2, textvariable=width_var, width=40).pack(side="left")

        row3 = ctk.CTkFrame(frame, fg_color="transparent")
        row3.pack(fill="x", padx=5, pady=(0,5))
        lbl_readout = ctk.CTkLabel(row3, text="Signal: 0.00", font=("Arial", 14, "bold"), text_color=color)
        lbl_readout.pack(side="left")
        
        ctk.CTkButton(row3, text="X", width=30, fg_color="#aa0000", hover_color="#ff0000",
                      command=lambda f=frame: self.remove_element(f, self.sensor_channels)).pack(side="right")
        
        self.sensor_channels.append({
            "frame": frame, "name_var": name_var, "center_var": center_var, "width_var": width_var, 
            "sensor_type_var": sensor_type_var, "lbl_readout": lbl_readout, "color": color
        })

    def remove_element(self, frame, target_list):
        frame.destroy()
        target_list[:] = [item for item in target_list if item["frame"] != frame]

    def get_source_spectra(self):
        source_type = self.source_type_var.get()
        if source_type == "Blackbody (Halogen)":
            try:
                temp = float(self.source_temp_var.get())
            except ValueError:
                temp = 3000
            return blackbody_spectrum(self.wl, temp)

        if source_type == "MTE6114W-WRC LED (1460nm)":
            return gaussian_bandpass(self.wl, 1460, 103, 1.0)

        if source_type == "HPIR104 Thermal Emitter":
            return bounded_blackbody_spectrum(self.wl, 903.15, min_nm=2000, max_nm=11000)

        return np.ones_like(self.wl)

    def get_sensor_spectra(self, sensor_type):
        if sensor_type == "InGaAs":
            return ingaas_responsivity(self.wl)
        if sensor_type == "InAsSb (2.7-5.3um)":
            return inassb_responsivity(self.wl)
        if sensor_type == "MCT (MIR)":
            return mct_responsivity(self.wl)
        return np.ones_like(self.wl)

    def get_channel_definitions(self):
        channels = []
        for ch in self.sensor_channels:
            try:
                center = float(ch["center_var"].get())
                width = float(ch["width_var"].get())
            except ValueError:
                continue

            channels.append({
                "name": ch["name_var"].get(),
                "center": center,
                "width": width,
                "sensor_type": ch["sensor_type_var"].get(),
            })

        return channels

    # --- Pop-up Viewer ---

    def show_material_spectra(self, mat_name):
        """Spawns a popup window showing the isolated spectral profile of a given material."""
        if mat_name not in self.material_library:
            return
            
        data = self.material_library[mat_name]
        
        # Create Popup
        popup = ctk.CTkToplevel(self)
        popup.title(f"Material Profile: {mat_name}")
        popup.geometry("600x450")
        popup.attributes("-topmost", True) # Force to top initially so it doesn't get lost
        
        # Setup Figure with Dual Axes
        fig, ax_alpha = plt.subplots(figsize=(6, 4), facecolor='#2b2b2b')
        ax_alpha.set_facecolor('#2b2b2b')
        ax_alpha.set_title(f"{mat_name} Base Spectral Properties", color='white')
        ax_alpha.set_xlabel("Wavelength (nm)", color='white')
        ax_alpha.tick_params(colors='white')
        
        # Plot Absorption on Left Axis
        color_alpha = '#ff3366'
        ax_alpha.set_ylabel("Absorption Coefficient (mm⁻¹)", color=color_alpha)
        ax_alpha.plot(self.wl, data["alpha"], color=color_alpha, label="Absorption (α)")
        ax_alpha.tick_params(axis='y', labelcolor=color_alpha)
        ax_alpha.grid(True, alpha=0.2)
        
        # Plot Refractive Index on Right Axis
        ax_n = ax_alpha.twinx()
        color_n = '#00ccff'
        ax_n.set_ylabel("Refractive Index (n)", color=color_n)
        
        # Handle constant vs array refractive index
        n_data = np.full_like(self.wl, data["n"]) if isinstance(data["n"], (int, float)) else data["n"]
        ax_n.plot(self.wl, n_data, color=color_n, linestyle='--', label="Refractive Index (n)")
        ax_n.tick_params(axis='y', labelcolor=color_n)
        
        # Combine Legends
        lines_1, labels_1 = ax_alpha.get_legend_handles_labels()
        lines_2, labels_2 = ax_n.get_legend_handles_labels()
        ax_alpha.legend(lines_1 + lines_2, labels_1 + labels_2, facecolor='#333333', edgecolor='white', labelcolor='white', loc='upper left')
        
        fig.tight_layout()
        
        # Embed in Tkinter
        canvas = FigureCanvasTkAgg(fig, master=popup)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

        marker_alpha, = ax_alpha.plot([], [], marker="o", color="white", markersize=6, linestyle="None")
        marker_n, = ax_n.plot([], [], marker="o", color="white", markersize=6, linestyle="None")
        annotation = ax_alpha.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            color="white",
            bbox=dict(boxstyle="round,pad=0.3", fc="#333333", ec="white", alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="white"),
        )
        annotation.set_visible(False)

        def on_click(event):
            if event.xdata is None:
                return

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
            else:
                return

            annotation.set_visible(True)
            canvas.draw_idle()

        canvas.mpl_connect("button_press_event", on_click)
        canvas.draw()

    def show_all_material_spectra(self):
        """Shows all material absorption curves in one popup."""
        popup = ctk.CTkToplevel(self)
        popup.title("Material Absorbance Comparison")
        popup.geometry("900x650")
        popup.attributes("-topmost", True)

        fig, (ax_raw, ax_norm) = plt.subplots(2, 1, figsize=(9, 6), facecolor='#2b2b2b', sharex=True)

        colors = {
            "PE": "#00ffcc",
            "EVOH": "#ffcc00",
            "Nylon": "#cc33ff",
            "Water": "#33ccff",
            "Air": "#888888",
        }

        for ax in (ax_raw, ax_norm):
            ax.set_facecolor('#2b2b2b')
            ax.tick_params(colors='white')
            ax.grid(True, alpha=0.2)
            for spine in ax.spines.values():
                spine.set_color('gray')

        plotted = []
        for mat_name, data in self.material_library.items():
            if mat_name == "Air":
                continue

            alpha = np.asarray(data["alpha"], dtype=float)
            color = colors.get(mat_name, None)
            raw_line, = ax_raw.plot(self.wl, alpha, label=mat_name, color=color, linewidth=1.8)

            max_alpha = np.nanmax(alpha)
            if max_alpha > 0:
                norm_alpha = alpha / max_alpha
                norm_line, = ax_norm.plot(self.wl, norm_alpha, label=mat_name, color=color, linewidth=1.8)
                plotted.append({
                    "name": mat_name,
                    "raw": alpha,
                    "norm": norm_alpha,
                    "raw_line": raw_line,
                    "norm_line": norm_line,
                })

        ax_raw.set_title("Absorption Coefficient", color='white')
        ax_raw.set_ylabel("alpha (mm^-1)", color='white')
        ax_raw.set_yscale("symlog", linthresh=0.01)

        ax_norm.set_title("Normalized Band Shapes", color='white')
        ax_norm.set_xlabel("Wavelength (nm)", color='white')
        ax_norm.set_ylabel("Relative alpha", color='white')

        for ax in (ax_raw, ax_norm):
            ax.legend(facecolor='#333333', edgecolor='white', labelcolor='white', loc='upper right')

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=popup)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

        raw_marker, = ax_raw.plot([], [], marker="o", color="white", markersize=6, linestyle="None")
        norm_marker, = ax_norm.plot([], [], marker="o", color="white", markersize=6, linestyle="None")
        raw_annotation = ax_raw.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            color="white",
            bbox=dict(boxstyle="round,pad=0.3", fc="#333333", ec="white", alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="white"),
        )
        norm_annotation = ax_norm.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            color="white",
            bbox=dict(boxstyle="round,pad=0.3", fc="#333333", ec="white", alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="white"),
        )
        raw_annotation.set_visible(False)
        norm_annotation.set_visible(False)

        def on_click(event):
            if event.xdata is None or event.inaxes not in (ax_raw, ax_norm):
                return

            idx = int(np.argmin(np.abs(self.wl - event.xdata)))
            x_val = self.wl[idx]
            series_key = "raw" if event.inaxes == ax_raw else "norm"
            click_xy = np.array([event.x, event.y])

            nearest = None
            nearest_distance = float("inf")
            for item in plotted:
                y_val = item[series_key][idx]
                pixel_xy = event.inaxes.transData.transform((x_val, y_val))
                distance = np.linalg.norm(pixel_xy - click_xy)
                if distance < nearest_distance:
                    nearest = item
                    nearest_distance = distance

            if nearest is None:
                return

            y_val = nearest[series_key][idx]
            if event.inaxes == ax_raw:
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

    def show_channel_matrix(self):
        """Shows effective absorption coefficients for the current channels."""
        channels = self.get_channel_definitions()
        materials = [name for name in self.material_library.keys() if name != "Air"]

        if not channels:
            return

        source_spectra = self.get_source_spectra()
        matrix = np.zeros((len(channels), len(materials)), dtype=float)
        channel_weights = []

        for row_idx, channel in enumerate(channels):
            filter_spectra = gaussian_bandpass(self.wl, channel["center"], channel["width"], 1.0)
            sensor_spectra = self.get_sensor_spectra(channel["sensor_type"])
            weight = source_spectra * filter_spectra * sensor_spectra
            weight_area = np.trapezoid(weight, self.wl)
            channel_weights.append(weight_area)

            for col_idx, material in enumerate(materials):
                alpha = np.asarray(self.material_library[material]["alpha"], dtype=float)
                if weight_area > 0:
                    matrix[row_idx, col_idx] = np.trapezoid(weight * alpha, self.wl) / weight_area
                else:
                    matrix[row_idx, col_idx] = 0.0

        singular_values = np.linalg.svd(matrix, compute_uv=False)
        nonzero_singular_values = singular_values[singular_values > 1e-12]
        rank = int(np.linalg.matrix_rank(matrix, tol=1e-9))
        if len(nonzero_singular_values) >= 2:
            condition = nonzero_singular_values[0] / nonzero_singular_values[-1]
        elif len(nonzero_singular_values) == 1:
            condition = float("inf")
        else:
            condition = float("inf")

        popup = ctk.CTkToplevel(self)
        popup.title("Channel Matrix")
        popup.geometry("1100x760")
        popup.attributes("-topmost", True)

        header = ctk.CTkLabel(
            popup,
            text="Effective alpha matrix: rows are sensor/filter channels, columns are materials",
            font=("Arial", 16, "bold"),
        )
        header.pack(pady=(10, 4))

        summary_text = (
            f"Channels: {len(channels)}    Materials: {len(materials)}    "
            f"Rank: {rank}    Condition: {condition:.3g}"
        )
        ctk.CTkLabel(popup, text=summary_text).pack(pady=(0, 8))

        fig, ax = plt.subplots(figsize=(8, 4.8), facecolor='#2b2b2b')
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
        ax.set_xticklabels(materials, color="white")
        ax.set_yticklabels([f"{ch['name']} ({ch['center']:.0f} nm)" for ch in channels], color="white")
        ax.tick_params(colors="white")
        for spine in ax.spines.values():
            spine.set_color("gray")

        for row_idx in range(len(channels)):
            for col_idx in range(len(materials)):
                value = matrix[row_idx, col_idx]
                ax.text(col_idx, row_idx, f"{value:.3g}", ha="center", va="center", color="white", fontsize=8)

        cbar = fig.colorbar(im, ax=ax)
        cbar.ax.tick_params(colors="white")
        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=popup)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=5)
        canvas.draw()

        table_lines = []
        table_lines.append("Effective alpha values are in mm^-1.")
        table_lines.append("Use absorbance A = -ln(signal/reference). Then approximately A = K * thickness.")
        table_lines.append("")
        table_lines.append("Channel".ljust(28) + "".join(material.rjust(14) for material in materials) + "    Weight")
        for row_idx, channel in enumerate(channels):
            name = f"{channel['name']} {channel['center']:.0f}nm"
            values = "".join(f"{matrix[row_idx, col_idx]:14.5g}" for col_idx in range(len(materials)))
            table_lines.append(name[:27].ljust(28) + values + f"    {channel_weights[row_idx]:.3g}")

        table_lines.append("")
        if len(channels) < len(materials):
            table_lines.append("Warning: fewer channels than materials. Thickness solving is underdetermined.")
        elif rank < len(materials):
            table_lines.append("Warning: matrix is rank deficient. Some materials look too similar in these channels.")
        elif condition > 100:
            table_lines.append("Warning: high condition number. Small signal noise may cause large thickness errors.")
        else:
            table_lines.append("Matrix looks reasonably separable for these modeled curves.")

        table_lines.append("Prefer columns that look different from each other; similar columns are hard to separate.")

        text_box = ctk.CTkTextbox(popup, height=170, font=("Consolas", 12))
        text_box.pack(fill="x", padx=10, pady=(5, 10))
        text_box.insert("1.0", "\n".join(table_lines))
        text_box.configure(state="disabled")

    # --- Core Physics Engine ---

    def run_live_simulation(self):
        # 1. Prepare the Source Object
        source_spectra = self.get_source_spectra()
        source_obj = LightSource("Source", self.wl, source_spectra)

        # 2. Prepare the Web Stack Objects
        web_stack_objs = []
        for i, layer_data in enumerate(self.web_layers):
            mat = layer_data["mat_var"].get()
            try: 
                d = float(layer_data["thick_var"].get())
            except ValueError: 
                d = 0.0
            
            alpha = self.material_library[mat]["alpha"]
            n_data = self.material_library[mat]["n"]
            
            layer_obj = MaterialLayer(
                name=f"Layer_{i}_{mat}", 
                thickness=d, 
                raw_wavelengths=self.wl, 
                alpha_data=alpha, 
                n_data=n_data
            )
            web_stack_objs.append(layer_obj)

        # Clear axes
        for ax in (self.ax_top, self.ax_bot):
            ax.clear()
            ax.set_facecolor('#2b2b2b')
            ax.tick_params(colors='white')
            for spine in ax.spines.values(): 
                spine.set_color('gray')

        self.ax_top.set_title("System Optics: Source Emittance & Web Transmission", color='white')
        self.ax_top.set_ylabel("Transmission / Intensity", color='white')
        self.ax_bot.set_title("Filtered Signals Reaching Detectors", color='white')
        self.ax_bot.set_xlabel("Wavelength (nm)", color='white')

        # Plot source output
        self.ax_top.plot(self.wl, source_spectra, color='#aaaaaa', linestyle="--", label="Source Output")

        first_channel = True

        for ch in self.sensor_channels:
            try:
                c_wl = float(ch["center_var"].get())
                fwhm = float(ch["width_var"].get())
            except ValueError:
                c_wl, fwhm = 2000, 10
            
            # 3. Prepare Filter and Sensor Objects
            filter_spectra = gaussian_bandpass(self.wl, c_wl, fwhm, 1.0)
            sensor_spectra = self.get_sensor_spectra(ch["sensor_type_var"].get())

            filter_obj = OpticalFilter(f"Filter_{c_wl}", self.wl, filter_spectra)
            sensor_obj = Sensor("Sensor", self.wl, sensor_spectra)

            # 4. RUN SIMULATION USING SHARED ENGINE
            results = run_simulation(self.wl, source_obj, web_stack_objs, filter_obj, sensor_obj)
            
            final_signal = results["final_signal"]
            channel_spectra = results["spectra"]["signal_spectrum"]
            
            # Plot the transmitted spectrum (web transmission) just once
            if first_channel:
                # Multiply bulk and interface transmission for the total realistic curve
                T_total = results["spectra"]["bulk_transmission"] * results["spectra"]["interface_transmission"]
                self.ax_top.plot(self.wl, T_total, color='white', label="Transmitted Spectrum (Bulk + Fresnel)")
                first_channel = False

            # Update UI Readout
            ch["lbl_readout"].configure(text=f"Signal: {final_signal:.2f}")

            name = ch["name_var"].get()
            color = ch["color"]
            
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
