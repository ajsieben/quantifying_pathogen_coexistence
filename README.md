# Quantifying the mechanisms of coexistence in disease ecology.

Repository of final code to accompany the manuscript "Quantifying the mechanisms of coexistence in disease ecology". The code is organized as follows:

**decomposition.R:** main script that calls the other model components and model functions, forming the overall structure of the model decompositions. The script calls parameter values from the "model_functions.R" script, including the time frame over which GRWR are averaged and decompositions are calculated.

**between_host_model.R:** code that runs pathogen transmission between hosts, host demography, and when within-host dynamics occur during each between-host time step.

**within_host_model.R** code that runs within-host infection dynamics following a series of ODEs modeling modified Lotka-Volterra dynamics.

**model_functions.R:** code with decomposition functions called in "decompositions.R" and other miscellaneous functions, roughly in order they are used in "decomposition.R". Parameter values for model runs are taken from "parameters.R" and used in "model_functions.R".

**parameters.R:** script with parameter values for different model simulations (e.g. neutral scenario, comp/col scenario, sensitivity analysis). Supplemental Figure 1 (varying within-host time-steps) based on neutral scenario parameters.

**analysis.R:** code used to generate figures and perform sensitivity analysis.

**rmpi.R:** script used for HPC at University of Wyoming.





