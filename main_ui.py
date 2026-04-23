import customtkinter as ctk
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import os 

# Import the backend components
from Components import LightSource, OpticalFilter, Sensor, MaterialLayer

# --- Physics & Data Helpers ---
def gaussian_peak(wl, center, width, amplitude):
    return amplitude * np.exp(-((wl - center) ** 2) / (2 * width ** 2))

def blackbody_spectrum(wl_nm, temp_k):
    h, c, k = 6.626e-34, 3.0e8, 1.38e-23
    wl_m = wl_nm * 1e-9
    intensity = (2 * h * c**2 / wl_m**5) / (np.exp(h * c / (wl_m * k * temp_k)) - 1)
    return intensity / np.max(intensity)

def ingaas_responsivity(wl):
    return np.maximum(np.interp(wl, [1500, 2000, 2300, 2500, 2600], [0.6, 0.9, 1.0, 0.8, 0.0]), 0)


import os
import numpy as np

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
        self.wl = np.arange(2000, 4000, 1.0)
        self.web_layers = []
        self.sensor_channels = []
        
        self.colors = ['#00ffcc', '#ff3366', '#ffcc00', '#cc33ff', '#33ccff']


        #   --- Pre-loaded Material Library ---
        self.material_library = {
            "Air": {"alpha": np.zeros_like(self.wl), "n": 1.0},
            "EVOH": {"alpha": gaussian_peak(self.wl, 1930, 30, 1.8) + gaussian_peak(self.wl, 2310, 25, 0.5), "n": 1.52},
            "Nylon": {"alpha": gaussian_peak(self.wl, 1530, 25, 1.2) + gaussian_peak(self.wl, 2050, 30, 1.5), "n": 1.53}
        }

        # Load both n and k from your new file
        # (Make sure to change the string to whatever you named the file!)
        pe_n, pe_k = load_database_file("PE_data.yml", self.wl)
        
        # Convert Extinction Coefficient (k) to Absorption (alpha) in mm^-1
        if pe_k is not None:
            # self.wl is in nanometers. Multiply by 1e-6 to get millimeters
            pe_alpha = (4 * np.pi * pe_k) / (self.wl * 1e-6)
        else:
            # Fallback to synthetic if k isn't found
            pe_alpha = gaussian_peak(self.wl, 1730, 20, 0.5) + gaussian_peak(self.wl, 2310, 25, 2.0)

        # Build the final PE material dictionary
        self.material_library["PE"] = {
            "alpha": pe_alpha,
            "n": pe_n if pe_n is not None else 1.51 
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
                          values=["Blackbody (Halogen)", "Flat Emission (Ideal)"]).pack(pady=5)
        
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
        ctk.CTkOptionMenu(row1, variable=sensor_type_var, values=["InGaAs", "Ideal (Flat)"], width=90).pack(side="right")

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
        canvas.draw()

    # --- Core Physics Engine ---

    def run_live_simulation(self):
        if self.source_type_var.get() == "Blackbody (Halogen)":
            try: temp = float(self.source_temp_var.get())
            except ValueError: temp = 3000
            source_spectra = blackbody_spectrum(self.wl, temp)
        else:
            source_spectra = np.ones_like(self.wl)
        
        T_bulk = np.ones_like(self.wl)
        for layer_data in self.web_layers:
            mat = layer_data["mat_var"].get()
            try: d = float(layer_data["thick_var"].get())
            except ValueError: d = 0.0
            
            alpha = self.material_library[mat]["alpha"]
            T_bulk *= np.exp(-alpha * d)

        base_intensity = source_spectra * T_bulk

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
        self.ax_top.plot(self.wl, T_bulk, color='white', label="Transmitted Spectrum")

        for ch in self.sensor_channels:
            try:
                c_wl = float(ch["center_var"].get())
                fwhm = float(ch["width_var"].get())
            except ValueError:
                c_wl, fwhm = 2000, 10
            
            filter_spectra = gaussian_peak(self.wl, c_wl, fwhm, 1.0)
            sensor_spectra = ingaas_responsivity(self.wl) if ch["sensor_type_var"].get() == "InGaAs" else np.ones_like(self.wl)

            channel_spectra = base_intensity * filter_spectra * sensor_spectra
            final_signal = np.trapz(channel_spectra, self.wl)
            
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
