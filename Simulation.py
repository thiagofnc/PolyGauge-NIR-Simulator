import numpy as np
from Components import LightSource, OpticalFilter, Sensor, MaterialLayer

def run_simulation(
    sim_wavelengths: np.ndarray,
    source: LightSource,
    web_stack: list[MaterialLayer],
    filter_obj: OpticalFilter,
    sensor: Sensor
) -> dict:
    """
    Simulates the entire optical path and calculates the final sensor signal.
    
    Parameters:
    - sim_wavelengths: The master array of wavelengths to run the simulation over.
    - source: LightSource object.
    - web_stack: A list of MaterialLayer objects representing the plastic film, in order.
    - filter_obj: OpticalFilter object.
    - sensor: Sensor object.
    
    Returns:
    - A dictionary containing the final integrated signal and intermediate spectral arrays.
    """
    
    # 1. Evaluate Source, Filter, and Sensor over the common wavelength axis
    I_0 = source.evaluate(sim_wavelengths)
    T_filter = filter_obj.evaluate(sim_wavelengths)
    S_sensor = sensor.evaluate(sim_wavelengths)

    # 2. Calculate Bulk Transmission (Beer-Lambert Law)
    T_bulk = np.ones_like(sim_wavelengths)
    for layer in web_stack:
        alpha = layer.alpha.evaluate(sim_wavelengths)
        # T = exp(-alpha * thickness)
        layer_transmission = np.exp(-alpha * layer.thickness)
        T_bulk *= layer_transmission

    # 3. Calculate Interface Transmission (Fresnel Equations)
    # We must account for the air on both sides of the web. Air n ≈ 1.0
    n_air = np.ones_like(sim_wavelengths)
    
    # Compile a list of all refractive indices in the path: Air -> Layers -> Air
    n_sequence = [n_air]
    for layer in web_stack:
        n_sequence.append(layer.n.evaluate(sim_wavelengths))
    n_sequence.append(n_air)

    T_interface = np.ones_like(sim_wavelengths)
    
    # Loop through every boundary between layers
    for i in range(len(n_sequence) - 1):
        n1 = n_sequence[i]
        n2 = n_sequence[i + 1]
        
        # Reflection coefficient for normal incidence
        R_boundary = ((n1 - n2) / (n1 + n2))**2
        # Transmission is 1 - Reflection
        T_interface *= (1 - R_boundary)

    # 4. Calculate Final Spectral Intensity reaching the sensor
    # Multiply all the continuous spectral curves together
    spectral_intensity = I_0 * T_interface * T_bulk * T_filter
    
    # 5. Calculate Final Sensor Response Spectrum
    sensor_response_spectrum = spectral_intensity * S_sensor

    # 6. Integrate to get the final measurable signal (e.g., total Amps or Volts)
    final_signal = np.trapezoid(sensor_response_spectrum, sim_wavelengths)

    return {
        "final_signal": final_signal,
        "spectra": {
            "wavelengths": sim_wavelengths,
            "source_intensity": I_0,
            "bulk_transmission": T_bulk,
            "interface_transmission": T_interface,
            "filter_transmission": T_filter,
            "signal_spectrum": sensor_response_spectrum
        }
    }