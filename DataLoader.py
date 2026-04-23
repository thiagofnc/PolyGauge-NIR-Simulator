import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

def load_real_material(filepath: str, sim_wavelengths_nm: np.ndarray) -> dict:
    """
    Loads n and k data from a standard optical database CSV and formats it 
    for the Web Gauging UI.
    
    Expected CSV format (RefractiveIndex.INFO standard):
    Wavelength (um), n, k
    1.50,            1.51, 0.0001
    ...
    """
    # 1. Load the raw data
    # Assuming CSV has headers: 'wl_um', 'n', 'k'
    df = pd.read_csv(filepath)
    
    raw_wl_nm = df['wl_um'].values * 1000.0  # Convert micrometers to nanometers
    raw_n = df['n'].values
    raw_k = df['k'].values

    # 2. Interpolate the raw data to match our UI's master wavelength array
    # We use fill_value="extrapolate" in case the DB data stops at 2500nm but our UI goes to 2600nm
    n_interp = interp1d(raw_wl_nm, raw_n, kind='linear', bounds_error=False, fill_value="extrapolate")
    k_interp = interp1d(raw_wl_nm, raw_k, kind='linear', bounds_error=False, fill_value="extrapolate")

    n_mapped = n_interp(sim_wavelengths_nm)
    k_mapped = k_interp(sim_wavelengths_nm)

    # 3. Convert Extinction Coefficient (k) to Absorption Coefficient (alpha)
    # Our UI uses thickness in mm, so we need alpha in mm^-1
    # wl_mm = wl_nm * 1e-6
    wl_mm = sim_wavelengths_nm * 1e-6
    
    # Avoid division by zero if wavelength array accidentally contains a 0
    with np.errstate(divide='ignore', invalid='ignore'):
        alpha_mapped = (4 * np.pi * k_mapped) / wl_mm

    # 4. Return the dictionary format that main_ui.py expects
    return {
        "alpha": alpha_mapped,
        "n": n_mapped
    }