# Micro risks and Pareto Improving policies 

The repository contains the code associated with the paper:

"Micro Risks and Pareto Improving Policies withLow Interest Rates" 
    
by Mark Aguiar, Manuel Amador and Cristina Arellano. 


## Structure

The subfolder `src` contains the main source code.

The subfolder `scripts` contains some of the analysis of the model for certain parameters. It contains both a julia script as well as the corresponding jupyter notebook. The scripts generate the figures and moments reported in the paper.  

The subfolder `output` contains the figures generates by the scripts, as well as some temporary calculations to speed up the simulations. 

## Running the code 

The code is in [Julia](https://julialang.org/downloads/).

To run the code, open a julia prompt at the root of this repository and type:

    julia> using Pkg 
    julia> Pkg.activate(".")
    julia> Pkg.instantiate()

The above will download the packages needed to run the code. 
  

To run the jupyter notebook, do:
  
    julia> using IJulia
    julia> notebook(dir=".")
  
That should open a browser with jupyter. Navigate to `scripts` to locate the notebooks. 

