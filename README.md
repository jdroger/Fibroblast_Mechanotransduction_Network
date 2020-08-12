# Fibroblast Chemo-/Mechano-transduction Model -- Development and Analysis Scripts

08.10.2020 Jesse Rogers

This repository provides analyses for a Netflux-generated ODE network of fibroblast signaling. Included here are scripts for performing model validation, mechano-chemo interaction, and network perturbation analyses, as well as current drug case studies and comprehensive drug screens.

## Required Programs/Toolboxes

- **MATLAB 2017a or newer**: needed to generate heatmap objects (MATLAB 2020a used for study)
  - *Parallel Computing Toolbox*: needed for drug screens on high performance computing environment
- **High performance computing environment/resource**: needed for comprehensive drug screens only
  - Provided batch scripts rely on Portable Batch System (PBS) job scheduler as implemented in the Palmetto cluster at Clemson University, and other environments may require alternative scripts

## Included Analyses

### 1. Model Validation

Qualitative validation of input-output and input-intermediate predictions with independent literature set.
- *ValidationSimulations.m*: Runs model simulations of network activation in response to single input stimuli, and categorizes changes in node activity compared to baseline conditions based on a 5% threshold
- *ValidationAnalysis.m*: Compares model-predicted changes to independent literature set and visualizes as a qualitative heatmap (Figure 2)

#### 1A. Parameter Sweep

Performs parameter sweep of EC50/n reaction parameters, predicting qualitative changes in node activity for each parameter set and comparing predictions to independent literature set to optimize parameters (Figure S2)
- *ValidationSimulations_paramsweep.m*: Repeats simulations for *ValidationSimulations.m* for all parameter sets
- *ValidationAnalysis_paramsweep.m*: Computes prediction accuracy for all parameter sets and visualizes top performing set

### 2. Mechano-Chemo Interaction Analysis

Simulation and comparison of dose-response behavior for biochemical inputs under varying levels of tension using changes in area-under-the-curve (dAUC) compared to baseline tension
- *InteractionSimulations.m*: Simulates steady-state network activation in response to 1% increments of each biochemical input under static levels of tension
- *InteractionAnalysis.m*: Calculates dAUC metrics from each dose-response curve following normalization step, and identifies reversal cases from AUC signs
- *InteractionVisualization.m*: Visualizes network-wide heatmaps of dAUC metrics (Figure 3A-B) and estimated density distributions for each input (Figure S3A)

### 3. Network Perturbation Analysis

Simulation of individual node knockdowns under low, medium, and high levels of tension, and identification of node sensitivity towards knockdown ("knockdown sensitivity") and influence on network-wide activity with knockdown ("knockdown influence")
- *PerturbationSimulations.m*: Simulates steady-state network activation in response to knockdown of individual nodes under a specified level of tension
  - Requires tension input reaction weight as argument (floating-point number between 0 and 1); levels 0.25, 0.5, and 0.75 used for study
- *PerturbationAnalysis.m*: Calculates knockdown sensitivity and influence from simulated activation levels, and identifies top 10 nodes for each metric
  - Requires tension input reaction weight as above
- *PerturbationVisualization.m*: Constructs heatmaps of top-ranking nodes in sensitivity and influence for each tension level (Figure 4A-C) and bar plots combining all tension levels (Figure 4D-E)

### 4. Drug Case Studies

Simulation of output/intermediate responses towards angiotensin receptor blocker (ARB) and neprilysin inhibitor (NEPi) drugs as measured in published experimental studies
- *DrugAnalysis.m*: Master script for running simulations of experimental studies, importing experimental data, processing simulation/experimental results, and generating all visualizations

### 5. Mechano-Adaptive Drug Screens

Comprehensive screening of all individual/combination drug targets for adaptive changes to matrix-related output expression
- *DrugScreen_individual.m*: Runs individual drug screen for both knockdowns (ko) and overexpression (oe)
- *DrugScreen_combinations_[type].m*: Runs combination drug screens for types dual ko [ko], dual overexpression [oe], node 1 ko + node 2 oe [ko-oe], and node 1 oe + node 2 ko [oe-ko]
- *DrugScreenAnalysis_individual.m*: Identifies mechano-adaptive perturbations from individual drug screen (Figure S6)
- *DrugScreenAnalysis_combinations.m*: Identifies mechano-adaptive perturbations from combinations drug screens (Figure 6)

#### 5A. HPC Batch Scripts
Includes bash shell scripts for submitting *DrugScreen* scripts to PBS scheduler. Recommended resources per script: 8 CPUs, 64GB memory, 48:00 runtime