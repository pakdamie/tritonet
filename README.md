This is the full code for the model and the subsequent analysis.



## The mathematical model
The full mathematical model is written in `rcpp` and should be in the code
`model_vectors_host.cpp`

## Simulation 
Any functions related to simulating dynamics from `model_vectors_host.cpp` is
in the simulate_functions.R.

## Calculation of the reproductive number
With the simulation output, these are then analyzed using functions
from the `simulate_functions.R` - notably for the reproductive
number.

## Plots
All main figure plots as well as the supplementary plots are listed in the 
`plotting_functions.R` - note that I use illustrator to deal with some very
minor labeling but nothing of the data is changed (I can easily use the unmodified
figures for the paper!)

