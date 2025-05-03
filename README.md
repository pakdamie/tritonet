This is the full code for the model simulation as well as the subsequent analysis.
This repository is by Damie Pak.

## The Main Infrastructure ("I want to get idea for how the code works!")

### The Mathematical Model
The full mathematical model is written in `rcpp` and should be in the code
`model_vectors_host.cpp`. If you're just interested in how I implement this
in Rcpp, I suggest looking just at this file. 

### Simulation 
Any functions related to simulating dynamics from `model_vectors_host.cpp` is
in the `simulate_functions.R`. You might be interested if you're trying to figure
out how I varied the different parameters.

### Calculation
With the simulation output, these are then analyzed using functions
from the `calculation_functions.R` - notably for the reproductive
number. Because I analytically derived the reproductive number, it's not
numerically simulated.

### Plots
All main figure plots as well as the supplementary plots are listed in the 
`plotting_functions.R` - note that I use illustrator to deal with some very
minor labeling but nothing of the data is changed (I can easily use the unmodified
figures for the paper!)

## Rebuilding and Reproducing  ("I want to reproduce your work!")

I provide a pipeline script R. It's called `main_pipeline.R`!
