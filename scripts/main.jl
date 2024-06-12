# -*- coding: utf-8 -*-
# # Micro Risks and Pareto Improving Policies 

@info "Initializing the environment..."
pkg_time = @elapsed begin 
    import Pkg; 
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end 
@info "Initializing the environment... done. Took $pkg_time seconds."


############################################################################
# # CLEANING UP THE OUTPUT DIRECTORY
############################################################################

@info "Cleaning up the output directory..."

# ## Delete all output files in the output directory 

for (root, dirs, files) in walkdir(joinpath(@__DIR__, "..", "output"))
    for file in files
        if endswith(file, ".pdf") || endswith(file, ".txt")
            rm(joinpath(root, file))
        end
    end
end

# ## Recreate the output directory structure if necessary 

mkpath(joinpath(@__DIR__, "..", "output", "tables"))
mkpath(joinpath(@__DIR__, "..", "output", "figures"))

@info "Cleaning up the output directory...done"

@info "Running the script..."
script_time = @elapsed include(joinpath(@__DIR__, "figures_and_tables.jl"))
@info "Running the script... done. Took $script_time seconds."

