# PolyGauge NIR Simulator

`PolyGauge NIR Simulator` is an early-stage Python project for modeling NIR/SWIR/MIR transmission through multilayer plastic films and exploring how narrowband filtered detector channels respond to that stack.

Right now, the repository contains two main ways to work with the model:

- A script-based simulation in `main.py` that builds a synthetic optical stack and runs a wavelength-by-wavelength calculation.
- A desktop GUI in `main_ui.py` that lets you interactively vary source settings, material layers, detector channels, and channel-matrix designs.

The project is clearly aimed at web gauging / film inspection workflows, especially polymer structures such as PE / EVOH / PE laminates or coextruded films.

## Current State

This is a solid prototype with a clear structure:

- Reusable spectral component classes are already separated into `Components.py`
- The core numerical pipeline is isolated in `Simulation.py`
- Data-loading logic exists in `DataLoader.py`
- A usable CustomTkinter front end is present in `main_ui.py`
- Real optical-constant data for polyethylene and water is bundled in the repo
- The GUI includes material curve comparison and channel-matrix tools for emitter/filter/sensor selection

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
- Can a set of channels mathematically separate PE, EVOH, Nylon 6, Nylon 66, and water thicknesses?

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
- material absorbance overlay viewer
- channel matrix viewer for checking whether the selected filters/sensors can separate material thicknesses

The GUI wavelength axis currently spans `1000-4000 nm`, which covers common NIR polymer bands plus the `3.0-3.5 um` stretch region. The detector choices are approximate:

- `InGaAs` for roughly `1500-2600 nm`
- `InAsSb (2.7-5.3um)` for MIR channels such as `3000` and `3400 nm`
- `MCT (MIR)` for roughly `2800-4000 nm`
- `Ideal (Flat)` for exploratory calculations without detector cutoff

Component-backed source options include:

- `MTE6114W-WRC LED (1460nm)`, modeled as a Gaussian source centered at `1460 nm` with `103 nm` spectral half-width from the Marktech 1460 nm emitter datasheet family.
- `HPIR104 Thermal Emitter`, modeled as a `903 K` near-blackbody emitter over its specified `2-11 um` output range.

For the filters in the current hardware shortlist, enter:

- `PIRF-INBP3400 Type 1`: center `3400 nm`, FWHM `140 nm`
- `Andover 3000.0 / 100.0 - 69350-B`: center `3000 nm`, FWHM `100 nm`
- No discrete filter with the 1460 nm LED: use the LED source option and a channel centered near `1460 nm`; a FWHM around `103 nm` represents the LED spectral width.

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

The GUI uses this file above `2600 nm`, so the strong PE `~3400 nm` C-H stretch region is based on measured optical constants rather than a hand-placed Gaussian.

### `PE_n_data.txt`

This looks like a simpler two-column refractive index table for PE:

- wavelength
- `n`

At the moment, the main GUI flow uses `PE_data.yml` directly, so this file appears to be auxiliary or legacy support data.

### `Water_data.yml`

This is liquid-water optical-constant data from the refractiveindex.info database, based on Hale and Querry's 25 C water `n,k` table from 0.2-200 micrometers.

The GUI converts `k` to absorption coefficient `alpha` in `mm^-1`, the same internal unit used by the rest of the simulator.

## Material Data Notes

The GUI material library is set up for NIR/SWIR/MIR absorbance-band comparison over `1000-4000 nm`.

- `Water` uses measured `n,k` data from `Water_data.yml` when that file is present. Its strong NIR bands are near `1450 nm` and `1940 nm`, so water layers should usually be modeled as very thin layers, often microns rather than tenths of a millimeter.
- `PE` uses engineering-scale Gaussian bands from `1000-2600 nm`, then switches to measured `n,k` data from `PE_data.yml` above `2600 nm`. That means the strong `~3400 nm` PE band is a real curve from optical constants, not just a guessed point.
- `EVOH`, `Nylon 6`, and `Nylon 66` currently use engineering-scale Gaussian band models. The band centers are sourced from NIR/IR spectroscopy literature, but the amplitudes are not calibrated extinction coefficients. Use them to compare band placement and channel sensitivity, not to infer absolute thickness.
- Channels above about `2600 nm` need the `MCT (MIR)` or `Ideal (Flat)` detector selection; the default `InGaAs` response intentionally cuts off in that region.

Current modeled NIR bands:

| Material | Main modeled bands |
| --- | --- |
| PE | `1210`, `1730`, `1764`, `2310`, `2350 nm`; measured optical constants above `2600 nm`, including the `~3400 nm` C-H stretch |
| EVOH | `1410`, `2012`, `2092`, `2310 nm`; modeled O-H stretch region near `3000-3300 nm` |
| Nylon 6 | `1485`, `1535`, `2040`, `2300`, `2355 nm`; modeled N-H stretch near `3030 nm` and CH2 stretch near `3410-3500 nm` |
| Nylon 66 | `1515`, `1565`, `2025`, `2075`, `2295`, `2350 nm`; modeled N-H stretch near `3060 nm` and CH2 stretch near `3415-3505 nm` |
| Water | measured `n,k`; prominent NIR absorption near `1450` and `1940 nm` |

Useful sources:

- Hale and Querry, "Optical constants of water in the 200-nm to 200-um wavelength region," Applied Optics 12, 555-563 (1973), via refractiveindex.info.
- Iwamoto et al., "FT-NIR Spectroscopic Study of OH Groups in Ethylene-Vinyl Alcohol Copolymer," Applied Spectroscopy 55, 864-870 (2001).
- Nylon 6 NIR studies in Journal of Molecular Structure report crystalline/amorphous Nylon 6 bands around `1485-1535 nm`, N-H combination bands near `~2040 nm`, and CH2 combination features around `2300-2355 nm`.
- Nylon 66 NIR studies, including work on heat-set carpet yarns and industrial polyamide monitoring, show related but shifted amide/N-H and methylene bands; the current Nylon 66 curve is a separate modeled polyamide response rather than a measured optical-constant dataset.
- PE NIR studies report CH2-related absorption around `1730-1764 nm` and strong hydrocarbon combination bands around `2300-2350 nm`.

## Viewing Curves

The GUI has two material-curve views:

- Click the eye button on an individual layer to view that material's absorption and refractive-index curves.
- Click `Compare Absorbance Curves` to overlay all material absorption curves.

Both popups support clicking on the plot. A marker appears at the nearest sampled wavelength and reports the wavelength plus the absorption value. The comparison popup includes:

- raw absorption coefficient curves in `mm^-1`
- normalized curves so weaker bands can still be compared by shape and wavelength

## Channel Matrix Workflow

Thickness solving should be done in absorbance space:

`A = -ln(signal / reference)`

For narrow enough channels, each channel approximately follows:

`A_channel = alpha_PE * d_PE + alpha_EVOH * d_EVOH + alpha_Nylon6 * d_Nylon6 + alpha_Nylon66 * d_Nylon66 + alpha_Water * d_Water`

The GUI's `CHANNEL MATRIX` button estimates the effective `alpha` value for each material in each current sensor/filter channel by weighting the material absorption curve by:

`source * filter * sensor`

Rows are channels, columns are materials, and values are effective `alpha` in `mm^-1`.

Use it to choose channels:

- You need at least as many independent channels as unknown material thicknesses.
- Similar columns mean two materials look alike to your selected channels and will be hard to separate.
- A high condition number means small detector noise can create large thickness errors.
- A zero or very low channel weight usually means the selected detector cannot see that wavelength range.

Practical workflow:

1. Add candidate sensor/filter channels in the left panel.
2. Use `InGaAs` for NIR channels and `InAsSb (2.7-5.3um)` or `MCT (MIR)` for `3000-3400 nm` channels.
3. Click `Compare Absorbance Curves` to inspect whether the chosen wavelengths sit on useful bands.
4. Click `CHANNEL MATRIX` to see whether the chosen channel set separates the materials.
5. Keep at least one reference or weak-absorption channel to normalize source/detector drift.

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
- MCT responsivity curve approximated with piecewise interpolation
- EVOH, Nylon 6, Nylon 66, and some PE NIR absorption bands modeled as Gaussians
- optical filter passbands modeled as Gaussian peaks

### Real / file-driven

- polyethylene `n` and `k` data from `PE_data.yml` above `2600 nm`
- liquid water `n` and `k` data from `Water_data.yml`

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
- opens material curve and channel-matrix popups from the center/right panel buttons

## Suggested Workflow

If you are developing this further, the most natural workflow is:

1. Use `Compare Absorbance Curves` to inspect candidate material bands.
2. Add candidate filter/sensor channels in `main_ui.py`.
3. Use `CHANNEL MATRIX` to check whether the selected channels can separate the unknown layer thicknesses.
4. Use `main.py` to validate backend math and inspect spectra when interface/Fresnel effects matter.
5. Move shared physics into `Simulation.py` so both script and GUI use the same engine.
6. Expand the material library with measured datasets for EVOH, Nylon 6, Nylon 66, and production PE grades.

## Known Limitations

This README should reflect the project honestly, so here are the main current limitations visible in the code:

- The GUI does not currently call the shared backend simulation function.
- Fresnel interface effects are modeled in `Simulation.py` but not in the GUI live simulation path.
- There is no dependency manifest file yet.
- There are no automated tests yet.
- There is no packaging, CLI, or installer yet.
- The material library is still small.
- EVOH, Nylon 6, and Nylon 66 absorption curves are modeled from band locations, not calibrated measured spectra.
- PE below `2600 nm` is still modeled; PE above `2600 nm` and water are file-driven measured optical constants.
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
5. Import measured EVOH, Nylon 6, and Nylon 66 absorbance or optical-constant curves from known-thickness samples.
6. Add export features for spectra, channel matrices, and channel results.
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
- bundled PE and water optical-constant datasets
- curve-comparison and channel-matrix tools for emitter/filter/sensor selection
- clear room to grow into a more serious engineering tool

For where it is right now, the architecture is already pointing in a good direction.
