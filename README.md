# Fibroblast Chemo-/Mechano-transduction Model -- Development and Analysis Scripts

08.10.2020 Jesse Rogers

This repository provides analyses for a Netflux-generated ODE network of fibroblast signaling. Included here are scripts for performing model validation, mechano-chemo interaction, and network perturbation analyses, as well as current drug case studies and comprehensive drug screens.

## Required Programs/Toolboxes

---

- **MATLAB 2017a or newer**: needed to generate heatmap objects
  - *Parallel Computing Toolbox*: needed for drug screens on high performance computing environment
- **High performance computing environment/resource**: needed for comprehensive drug screens only
  - Provided batch scripts rely on Portable Batch System (PBS) job scheduler as implemented in the Palmetto cluster at Clemson University, and other resources may require alternative scripts

## Included Analyses

---

- **1_ModelValidation**: Qualitative validation of input-output and input-intermediate predictions with independent literature set. Consists of functions for model simulations (validationsimulations.m) and comparison with literature (validationanalysis.m). Generates Figure 2
  - *1A_ParameterSweep*: Sweep of EC50/n parameters for optimizing prediction accuracy. Generates Supplemental Figure S2