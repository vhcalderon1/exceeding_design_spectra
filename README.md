<!--
SPDX-FileCopyrightText: 2025 Stanford University

SPDX-License-Identifier: MIT
-->

# Observed ground motions that exceeded design response spectra in the Western United States

This repository contains MATLAB code, data-processing scripts, and helper functions accompanying:

> Victor H. Calderon & Jack W. Baker (2025).  
> "Observed ground motions that exceeded design response spectra in the Western United States".  
> *Submitted to Earthquake Spectra* (in review).

All code is designed for reproducible, publication-quality seismic ground motion analysis. The project identifies and analyzes ground motions in the Western U.S. that exceeded modern seismic design response spectra, as described in the above paper.

---

## Overview

The scripts in this repository perform the following key tasks:
- **Screen set of NGA-West2 recordings** to find ground motions that exceed DE, MCE_R, probablisitic MCE_R and UHS target spectra.
- **Generate tables and figures** quantifying the number, magnitude, distance, and period ranges of exceedances.
- **Compare observed spectra to target and disaggregated values** for physical insight.

Sample outputs are written to `/results/figures/`, ready for direct inclusion in publications.

---

## Repository Structure

