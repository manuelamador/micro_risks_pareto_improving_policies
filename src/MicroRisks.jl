module MicroRisks

using LinearAlgebra
using ProgressMeter
using Roots
using NLsolve
using StructArrays
using StatsBase
using Tullio
using PrettyTables

include("types.jl")
include("labor_preferences.jl")
include("technology_and_fiscal.jl")
include("calibration.jl")
include("auxiliary_functions.jl")
include("household_problem_and_aggregation.jl")
include("stationary_equilibria.jl")
include("transition.jl")
include("statistics.jl")
include("jacobian.jl")

const _TOL = 1e-12
const _ZERO_FTOL = 1e-08
const _MAX_ITERS = 20_000
const _ZERO_MAX_ITERS = 300

export 

# types.jl 
    CRRA, EZ, GHH, FixedLabor, Household, HouseholdWorkspace, CobbDouglasTechnology, JacobianCache,

# calibration.jl
    calibration, 

# technology.jl 
    golden_rule_k, output, y, 

# household_problem_and_aggregation
    stationary, stationary!, is_pol_valid, 
    aggregate_labor, aggregate_assets, consumption_alloc, aggregate_consumption, 

# stationary_equilibria.jl
    stationary_laissez_faire, stationary_equilibrium_given_k_b,

# transition.jl
    solve_RPI_transition,

# statistics_and_plots.jl
    generate_tables, 

# jacobian.jl 
    pv_elasticities, pv_elasticities!, jacobian_row, jacobian,

# auxiliary_functions.jl
    grid

end