import numpy as np
import matplotlib.pyplot as plt
from Components import LightSource, OpticalFilter, Sensor, MaterialLayer
from Simulation import run_simulation

# --- Helper Functions for Synthetic Data ---

def gaussian_peak(wavelengths, center, width, amplitude):
    """Generates a synthetic absorption peak."""
    return amplitude * np.exp(-((wavelengths - center) ** 2) / (2 * width ** 2))

def blackbody_spectrum(wavelengths_nm, temp_k=3000):
    """Approximates a Halogen lamp using Planck's Law."""
    h = 6.626e-34
    c = 3.0e8
    k = 1.38e-23
    wl_m = wavelengths_nm * 1e-9
    intensity = (2 * h * c**2 / wl_m**5) / (np.exp(h * c / (wl_m * k * temp_k)) - 1)
    # Normalize to a peak of 1 for easier visualization
    return intensity / np.max(intensity)

def ingaas_responsivity(wavelengths_nm):
    """Approximates an extended InGaAs sensor response (Amps/Watt)."""
    # Ramps up, peaks around 2100nm, cuts off sharply around 2550nm
    response = np.interp(wavelengths_nm, 
                         [1500, 2000, 2300, 2500, 2600], 
                         [0.6,  0.9,  1.0,  0.8,  0.0])
    return np.maximum(response, 0)

# --- 1. Set Up the Master Wavelength Array ---
# 1500 nm to 2600 nm in 1 nm increments (Standard NIR gauging range)
wl = np.arange(1000, 4000, 1.0)

# --- 2. Generate Synthetic Test Data ---

# Light Source (Tungsten Halogen)
source_data = blackbody_spectrum(wl)
source = LightSource("Halogen_Lamp", wl, source_data)

# Sensor (InGaAs)
sensor_data = ingaas_responsivity(wl)
sensor = Sensor("InGaAs_Detector", wl, sensor_data)

# Optical Filters
# A "Measure" filter targeted directly at the PE C-H peak (2310 nm)
# A "Reference" filter targeted at a valley where absorption is low (2200 nm)
measure_filter_data = gaussian_peak(wl, center=2310, width=15, amplitude=0.85) # 85% peak transmission
ref_filter_data = gaussian_peak(wl, center=2200, width=15, amplitude=0.85)
measure_filter = OpticalFilter("PE_Measure_Filter", wl, measure_filter_data)
reference_filter = OpticalFilter("Reference_Filter", wl, ref_filter_data)

# Material Absorption Coefficients (units: mm^-1)
# PE: Peaks at 1730nm and 2310nm
alpha_pe = gaussian_peak(wl, 1730, 20, 0.5) + gaussian_peak(wl, 2310, 25, 2.0)
# EVOH: Peak at 1930nm (OH), plus a smaller CH peak
alpha_evoh = gaussian_peak(wl, 1930, 30, 1.8) + gaussian_peak(wl, 2310, 25, 0.5)
# Nylon: Peaks at 1530nm and 2050nm (NH)
alpha_nylon = gaussian_peak(wl, 1530, 25, 1.2) + gaussian_peak(wl, 2050, 30, 1.5)

# --- 3. Build the Optical Stack ---
# Remember: If alpha is in mm^-1, thickness must be in mm!
# Let's simulate a 5-layer coextruded film: PE / EVOH / PE (ignoring tie layers for simplicity)
layer1_pe = MaterialLayer("PE_Outer", thickness=0.050, raw_wavelengths=wl, alpha_data=alpha_pe, n_data=1.51)
layer2_evoh = MaterialLayer("EVOH_Core", thickness=0.015, raw_wavelengths=wl, alpha_data=alpha_evoh, n_data=1.52)
layer3_pe = MaterialLayer("PE_Inner", thickness=0.050, raw_wavelengths=wl, alpha_data=alpha_pe, n_data=1.51)

web_stack = [layer1_pe, layer2_evoh, layer3_pe]

# --- 4. Run the Simulations ---

# Run once with the Measure Filter
results_measure = run_simulation(wl, source, web_stack, measure_filter, sensor)

# Run once with the Reference Filter
results_reference = run_simulation(wl, source, web_stack, reference_filter, sensor)

# --- 5. Output and Visualization ---

print("=== Simulation Results ===")
print(f"Signal through Measure Filter (2310nm): {results_measure['final_signal']:.4f} units")
print(f"Signal through Reference Filter (2200nm): {results_reference['final_signal']:.4f} units")

# The ratio is what web gauging systems actually use to determine thickness!
ratio = results_measure['final_signal'] / results_reference['final_signal']
print(f"Signal Ratio (Measure/Reference): {ratio:.4f}")

# Plotting the physics
plt.figure(figsize=(12, 8))

# Subplot 1: Stack Transmission
plt.subplot(2, 1, 1)
plt.plot(wl, results_measure['spectra']['bulk_transmission'], label='Bulk Transmission (Beer-Lambert)', color='black')
plt.plot(wl, results_measure['spectra']['interface_transmission'], label='Interface Transmission (Fresnel)', color='gray', linestyle='--')
plt.title("Web Spectral Transmission (PE/EVOH/PE)")
plt.ylabel("Transmittance (0 to 1)")
plt.legend()
plt.grid(True, alpha=0.3)

# Subplot 2: The Filtered Signals hitting the Sensor
plt.subplot(2, 1, 2)
# We plot the signal spectrum *before* integrating it to a single number
plt.plot(wl, results_measure['spectra']['signal_spectrum'], label='Detector Signal (Measure Filter)', color='blue')
plt.plot(wl, results_reference['spectra']['signal_spectrum'], label='Detector Signal (Reference Filter)', color='orange')
plt.title("Final Spectral Signal on Detector")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Signal Intensity (Relative)")
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()