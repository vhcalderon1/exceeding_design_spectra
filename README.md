<!--
SPDX-FileCopyrightText: 2025 Stanford University

SPDX-License-Identifier: MIT
-->

# Observed ground motions that exceeded design response spectra in the Western United States

The project identifies and analyzes ground motions in the Western U.S. that exceeded modern seismic design response spectra, as described in:

> Victor H. Calderon & Jack W. Baker (2025).  
> "Observed ground motions that exceeded design response spectra in the Western United States".  
> *Submitted to Earthquake Spectra* (in review).

This repository contains MATLAB code, data-processing scripts, and helper functions.
If you encountered any issue, feel free to contact vhcalderon (at) princeton (dot) edu

## Overview

The scripts in this repository perform the following key tasks:
- **Screen set of NGA-West2 recordings** to find ground motions that exceed DE, MCE_R, probablisitic MCE_R and UHS target spectra.
- **Generate tables and figures** quantifying the number, magnitude, distance, and period ranges of exceedances.
- **Compare observed spectra to target and seismic hazard disaggregated values** for physical insight.

Sample outputs are written to `/results/figures/`, ready for direct inclusion in publications.

## Repository Structure
```
exceeding_design_spectra/
│
├─ main.m                 % reproducible entry-point (see below)
├─ scripts/               % high-level analyses & plotting
├─ src/                   % helper functions (spectra, interpolation, etc.)
├─ data/
│   ├─ raw/               % large external files (not version-controlled)
│   └─ (derived)          % .mat & .xlsx outputs created by scripts
└─ results/
    └─ figures/           % publication-ready PDFs/PNGs
```
## Getting Started

Clone this repository and install any necessary dependencies (see comments in scripts).

Most analysis can be performed by running the `main.m` script from the project root. Make sure you have all required data files in the `data/` folder. Figures and results will be saved to the `results/` folder.

## Key Scripts

The following scripts are run automatically by `main.m`:

- `disagg_residual_UHS.m`  
  Disaggregation of UHS residuals (see Fig. 9 of the paper)

- `disagg_mag_dist_UHS.m`  
  Magnitude–distance analysis for UHS exceedances (Fig. 9)
  
- `ratiosa_vs_tmax.m`  
  SA/UHS ratio vs. T/Tmax (Fig. 8)

- `magnitude_distance.m`  
  Magnitude–distance threshold plots (Fig. 5)

- `number_exceedances.m`  
  Number of spectral exceedances as a function of period (Fig. 4)

- `set_rsns_plot.m`  
  Map and plot of RSNs with exceedances (Fig. 2, 3, and 6)

- `usable_gms.m`  
  Usable ground motions vs. period (Fig. 1)

Additional batch, table, or data-preparation scripts are available in `scripts/` and can be run independently.

## Data

Large raw datasets must be downloaded manually. Specifically, gridded risk-targeted and 84th-percentile spectra must be obtained from [USGS ScienceBase](https://doi.org/10.5066/P9I0R4O6) and placed in `data/raw/Risk Target MCE Conterminous`.

## Licensing and Citation

All code in this repository is licensed under the [MIT License](LICENSE).

If you use these scripts or analysis in your research, please cite:

```bibtex
@article{calderonbaker2025observed,
  title   = {Observed ground motions that exceeded design response spectra in the Western United States},
  author  = {Victor H. Calderon and Jack W. Baker},
  year    = {2025},
  journal = {Earthquake Spectra (in review)},
  note    = {GitHub: https://github.com/vhcalderon1/exceeding_design_spectra}
}
```


