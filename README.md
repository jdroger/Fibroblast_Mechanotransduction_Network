# Fibroblast Chemo-/Mechano-transduction Model -- Development and Analysis Scripts

08.10.2020 Jesse Rogers

This repository provides analyses for a Netflux-generated ODE network of fibroblast signaling. Included here are scripts for performing model validation, mechano-chemo interaction, and network perturbation analyses, as well as current drug case studies and comprehensive drug screens.

## Required Programs/Toolboxes

- **MATLAB 2017a or newer**: needed to generate heatmap objects
  - *Parallel Computing Toolbox*: needed for drug screens on high performance computing environment
- **High performance computing environment/resource**: needed for comprehensive drug screens only
  - Provided batch scripts rely on Portable Batch System (PBS) job scheduler as implemented in the Palmetto cluster at Clemson University, and other environments may require alternative scripts

## Included Analyses

- **1_ModelValidation**: Qualitative validation of input-output and input-intermediate predictions with independent literature set.
  - *ValidationSimulations.m*: Runs model simulations of network activation in response to single input stimuli, and categorizes changes in node activity compared to baseline conditions based on a 5% threshold
  - *ValidationAnalysis.m*: Compares model-predicted changes to independent literature set and visualizes as a qualitative heatmap (Figure 2)
  - *1A_ParameterSweep*: Performs parameter sweep of EC50/n reaction parameters, predicting qualitative changes in node activity for each parameter set and comparing predictions to independent literature set to optimize parameters (Figure S2)
    - *ValidationSimulations_paramsweep.m*: Repeats simulations for *ValidationSimulations.m* for all parameter sets
    - *ValidationAnalysis_paramsweep.m*: Computes prediction accuracy for all parameter sets and visualizes top performing set
  