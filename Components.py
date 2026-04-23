import numpy as np
from scipy.interpolate import interp1d

class SpectralProfile:
    """
    Base class representing any property that varies across wavelengths.
    Handles the interpolation so all components can sync to a master wavelength array.
    """
    def __init__(self, name: str, raw_wavelengths: np.ndarray, data: np.ndarray):
        self.name = name
        # We use extrapolate so the simulation doesn't crash if arrays are slightly mismatched
        self._interp = interp1d(raw_wavelengths, data, kind='linear', 
                                bounds_error=False, fill_value="extrapolate")

    def evaluate(self, sim_wavelengths: np.ndarray) -> np.ndarray:
        """Returns the data array mapped to the simulation's common wavelength axis."""
        return self._interp(sim_wavelengths)


class LightSource(SpectralProfile):
    """Represents the IR emitter. Data is spectral irradiance (e.g., W/m^2/nm)."""
    pass


class OpticalFilter(SpectralProfile):
    """Represents a narrow-band filter. Data is Transmittance (0.0 to 1.0)."""
    pass


class Sensor(SpectralProfile):
    """Represents the photodetector. Data is spectral responsivity (e.g., Amps/Watt)."""
    pass


class MaterialLayer:
    """
    Represents a physical layer in the plastic web. 
    It doesn't inherit from SpectralProfile directly because a material requires 
    TWO spectral profiles (absorption and refractive index), plus a physical thickness.
    """
    def __init__(self, name: str, thickness: float, 
                 raw_wavelengths: np.ndarray, alpha_data: np.ndarray, n_data):
        
        self.name = name
        self.thickness = thickness
        
        # Profile for absorption coefficient (alpha)
        self.alpha = SpectralProfile(f"{name}_alpha", raw_wavelengths, alpha_data)
        
        # Profile for refractive index (n)
        if isinstance(n_data, (int, float)):
            # If n is constant, create a dummy profile that always returns that constant
            self.n = SpectralProfile(f"{name}_n", raw_wavelengths, np.full_like(raw_wavelengths, n_data, dtype=float))
        else:
            self.n = SpectralProfile(f"{name}_n", raw_wavelengths, n_data)