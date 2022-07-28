# Micro Risks and Pareto Improving Policies 

The repository contains the code associated with the paper:

["Micro Risks and (Robust) Pareto Improving Policies with Low Interest Rates"](https://manuelamador.me/files/rlessthang.pdf)
    
by Mark Aguiar, Manuel Amador and Cristina Arellano (2021).


## Structure

The subfolder `src` contains the main source code.

The subfolder `scripts` contains contains both a julia script as well as the corresponding jupyter notebook that generates the figures and moments reported in the paper.  

The subfolder `output` contains the figures generates by the script. 

## Running the code 

The code is written in [Julia](https://julialang.org/downloads/). The code uses multithreading if available, so make sure to start julia with the ability to run multiple threads. See [Julia Multi-threading](https://docs.julialang.org/en/v1/manual/multi-threading/#man-multithreading).

To run the code, open a julia prompt at the root of this repository. 
Then type:

    julia> using Pkg 
    julia> Pkg.activate(".")
    julia> Pkg.instantiate()

The above will download the Julia packages needed to run the programs.

To run the main script that generates the figures and statistics, do:
  
    julia> include(joinpath("scripts", "paper_figures_tables.jl")) 
    
The script takes about 50 minutes to run a in a server with AMD EPYC 7542 32-Core Processor with 20 cores available and 32GB of RAM.  

Associated with the script, there is a Jupyter notebook that can be run instead. Run [Jupyter](https://jupyter.org/) and navigate to `scripts` to locate the notebook and run it.

