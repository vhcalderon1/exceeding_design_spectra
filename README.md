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
├── main.m # Main driver script for end-to-end analysis
├── scripts/ # Top-level analysis and plotting scripts (run by main.m)
├── src/ # Helper functions (e.g. spectral interpolation, disaggregation, utility)
├── data/
│ ├── raw/ # External datasets (NGA-West2, USGS grids, not version controlled)
│ └── derived/ # Output .mat/.xlsx files from analysis scripts
└── results/
└── figures/ # Publication-quality plots (PDF, PNG)
```
