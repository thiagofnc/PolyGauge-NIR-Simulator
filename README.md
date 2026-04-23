# PolyGauge NIR Simulator

`PolyGauge NIR Simulator` is an early-stage Python project for modeling near-infrared transmission through multilayer plastic films and exploring how narrowband filtered detector channels respond to that stack.

Right now, the repository contains two main ways to work with the model:

- A script-based simulation in `main.py` that builds a synthetic optical stack and runs a wavelength-by-wavelength calculation.
- A desktop GUI in `main_ui.py` that lets you interactively vary source settings, material layers, and detector channels.

The project is clearly aimed at web gauging / film inspection workflows, especially polymer structures such as PE / EVOH / PE laminates or coextruded films.

## Current State

This is a solid prototype with a clear structure:

- Reusable spectral component classes are already separated into `Components.py`
- The core numerical pipeline is isolated in `Simulation.py`
- Data-loading logic exists in `DataLoader.py`
- A usable CustomTkinter front end is present in `main_ui.py`
- Real material data for polyethylene is bundled in the repo

It is not yet a packaged library or production application, but it already demonstrates the core physics and user workflow.

## What The Project Does

At a high level, the simulator models this signal chain:

1. A light source emits radiation across a wavelength range.
2. That light passes through a multilayer plastic web.
3. Each layer attenuates the spectrum according to its absorption coefficient and thickness.
4. Interfaces between materials reduce transmission due to Fresnel reflection.
5. A narrowband optical filter selects a region of the spectrum.
6. A detector weights the transmitted light by its responsivity.
7. The weighted spectrum is integrated to produce a final scalar signal.

That makes the project useful for questions like:

- How sensitive is a measurement channel near a polymer absorption band?
- How does a reference channel compare against a measurement channel?
- How does film stack design affect transmission?
- How do source spectrum, filter placement, and detector response interact?

## Repository Contents

### `main.py`

This is the scripted example and the clearest reference implementation of the backend flow.

It:

- Builds a master wavelength grid from `1000` to `3999 nm`
- Synthesizes a halogen-like source spectrum with a blackbody approximation
- Synthesizes an InGaAs detector responsivity curve
- Creates two Gaussian optical filters:
  - a PE measurement channel near `2310 nm`
  - a reference channel near `2200 nm`
- Builds a three-layer polymer stack:
  - `PE_Outer`
  - `EVOH_Core`
  - `PE_Inner`
- Runs the simulation for both filters
- Prints final scalar channel signals and their ratio
- Plots:
  - bulk transmission
  - interface transmission
  - detector signal spectra for both channels

This file is the best place to start if you want to understand the intended physical model.

### `Components.py`

This module defines the reusable optical objects used by the simulator.

#### `SpectralProfile`

Base class for any wavelength-dependent quantity. It stores raw spectral data and interpolates it onto a common simulation wavelength axis.

Used for:

- source intensity
- filter transmission
- detector responsivity
- material absorption
- material refractive index

#### `LightSource`

Represents the emitter spectrum.

#### `OpticalFilter`

Represents a bandpass or narrowband transmission curve.

#### `Sensor`

Represents detector responsivity versus wavelength.

#### `MaterialLayer`

Represents a physical film layer with:

- a name
- a thickness
- an absorption profile `alpha`
- a refractive-index profile `n`

This class is especially important because it bridges spectral data and physical thickness.

### `Simulation.py`

This is the backend physics engine currently used by `main.py`.

The `run_simulation(...)` function:

- evaluates all wavelength-dependent inputs on a common axis
- computes Beer-Lambert bulk attenuation for each layer
- computes Fresnel transmission losses at each interface
- multiplies source, web, filter, and sensor spectra together
- integrates the final detector response with `np.trapezoid`

It returns both:

- `final_signal`: a scalar measurement value
- `spectra`: intermediate arrays for debugging, plotting, and analysis

This separation is a strong design choice because it makes the model easier to inspect and reuse.

### `DataLoader.py`

This module loads real optical material data from CSV-style sources and converts extinction coefficient `k` into absorption coefficient `alpha`.

Specifically, it:

- reads `wl_um`, `n`, and `k`
- converts wavelength from micrometers to nanometers
- interpolates `n` and `k` to the simulation wavelength axis
- converts `k` to `alpha` using:

`alpha = (4 * pi * k) / wavelength`

with wavelength converted to millimeters so the resulting `alpha` is in `mm^-1`.

This is the bridge from external optical databases into the simulator’s internal material format.

### `main_ui.py`

This is the interactive desktop application built with `customtkinter` and Matplotlib.

The GUI currently provides:

- source selection
  - blackbody-like halogen source
  - flat ideal source
- source temperature entry
- editable sensor/filter channels
  - channel name
  - center wavelength
  - width
  - detector type
- editable web stack
  - selectable material
  - thickness per layer
- live spectral plots
- per-channel signal readouts
- popup material spectral viewer

Default startup configuration is already meaningful:

- source: blackbody halogen
- channels:
  - `Measure (PE)` around `2310 nm`
  - `Reference` around `2200 nm`
- layers:
  - `PE`
  - `EVOH`
  - `PE`

One important implementation detail:

The GUI currently contains its own simplified physics path instead of calling `run_simulation(...)` from `Simulation.py`.

That means:

- `main.py` includes both Beer-Lambert and Fresnel interface effects
- `main_ui.py` currently models bulk absorption only in its live update path

So the GUI is best viewed as an interactive prototype, while `main.py` is the more complete backend example.

### `PE_data.yml`

This is a bundled polyethylene optical dataset in a refractiveindex.info-style format.

From the file contents, it includes tabulated:

- wavelength
- refractive index `n`
- extinction coefficient `k`

The data appears to be used by the GUI to create a more realistic PE material instead of relying only on synthetic Gaussian peaks.

### `PE_n_data.txt`

This looks like a simpler two-column refractive index table for PE:

- wavelength
- `n`

At the moment, the main GUI flow is using `PE_data.yml` directly, so this file appears to be auxiliary or legacy support data.

## Physics Model Included So Far

The project already covers several useful pieces of optical modeling:

- wavelength-domain simulation
- interpolated spectral profiles
- Beer-Lambert absorption through stacked layers
- Fresnel losses at normal-incidence interfaces
- blackbody-like source approximation
- narrowband Gaussian filter approximation
- detector responsivity weighting
- integrated scalar detector output

This is enough to make the simulator educational and directionally useful for concept validation.

## What Is Synthetic vs Real

The repo currently mixes synthetic and real inputs.

### Synthetic / approximated

- halogen source modeled with Planck’s law approximation
- InGaAs responsivity curve approximated with piecewise interpolation
- several material absorption bands modeled as Gaussians
- optical filter passbands modeled as Gaussian peaks

### Real / file-driven

- polyethylene `n` and `k` data from `PE_data.yml`

That hybrid setup makes sense for a prototype: you get realistic structure where you have data, while keeping the rest fast to modify.

## Installation

There is no `requirements.txt` or `pyproject.toml` yet, so dependencies need to be installed manually.

Based on the imports in this repository, you currently need:

- `numpy`
- `matplotlib`
- `scipy`
- `pandas`
- `customtkinter`

Example:

```powershell
pip install numpy matplotlib scipy pandas customtkinter
```

## How To Run

### Run the scripted simulation

```powershell
python main.py
```

Expected behavior:

- prints measurement-channel signal
- prints reference-channel signal
- prints their ratio
- opens Matplotlib plots showing transmission and detector spectra

### Run the GUI

```powershell
python main_ui.py
```

Expected behavior:

- opens a desktop window
- lets you add/remove layers
- lets you add/remove channels
- recalculates channel signals
- updates plots of source/web transmission and filtered detector response

## Suggested Workflow

If you are developing this further, the most natural workflow is:

1. Use `main.py` to validate the backend math and inspect spectra.
2. Use `main_ui.py` to experiment with stack design and channel placement.
3. Move shared physics into `Simulation.py` so both script and GUI use the same engine.
4. Expand the material library with more real datasets.

## Known Limitations

This README should reflect the project honestly, so here are the main current limitations visible in the code:

- The GUI does not currently call the shared backend simulation function.
- Fresnel interface effects are modeled in `Simulation.py` but not in the GUI live simulation path.
- There is no dependency manifest file yet.
- There are no automated tests yet.
- There is no packaging, CLI, or installer yet.
- The material library is still small.
- Some data parsing in the GUI is custom and separate from `DataLoader.py`, so data-loading logic is duplicated.
- The project currently assumes relatively simple normal-incidence transmission and does not include scattering, angular effects, or instrument noise modeling.

## Recommended Next Steps

If you keep building this project, the highest-value improvements would likely be:

1. Refactor `main_ui.py` to use `Simulation.run_simulation(...)`.
2. Centralize all file loading through `DataLoader.py`.
3. Add a `requirements.txt` or `pyproject.toml`.
4. Add a few test cases for:
   - Beer-Lambert attenuation
   - Fresnel boundaries
   - material interpolation
   - signal integration
5. Expand the material library with real optical constants for more polymers and barrier layers.
6. Add export features for spectra and channel results.
7. Add calibration / regression tools to estimate thickness or composition from channel ratios.

## Why The New Name

The old visible title in the GUI was `NIR Web Gauging Simulator - Pro Edition`.

`PolyGauge NIR Simulator` is a better fit for the code currently in this repository because it is:

- shorter
- more specific to polymer film gauging
- easier to reuse as the project grows
- descriptive without sounding like a finished commercial product

## Summary

So far, this project is a promising NIR film-transmission and web-gauging prototype with:

- a reusable simulation backend
- an interactive GUI
- a bundled PE dataset
- clear room to grow into a more serious engineering tool

For where it is right now, the architecture is already pointing in a good direction.
