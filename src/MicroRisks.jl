module MicroRisks

using LinearAlgebra
using ProgressMeter
using Roots
using NLsolve
using StructArrays
using StatsBase
using Tullio
using PrettyTables
using Interpolations

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

#############################################################################################
# Aggregate risk types/methods
#

include(joinpath("aggregate_risk_extension", "types_aggregate_risk.jl"))
include(joinpath("aggregate_risk_extension","household_aggregate_risk.jl"))
include(joinpath("aggregate_risk_extension","transition_aggregate_risk.jl"))
include(joinpath("aggregate_risk_extension","technology_and_fiscal_aggregate_risk.jl"))

# CONSTANTS 

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
    stationary_laissez_faire, stationary_equilibrium_given_k_b, stationary_high_b,

# transition.jl
    solve_RPI_transition,

# statistics_and_plots.jl
    generate_tables, generate_tables_all, 

# jacobian.jl 
    pv_elasticities, pv_elasticities!, jacobian_row, jacobian,

# auxiliary_functions.jl
    grid,

#############################################################################################
# Aggregate risk components 
#

# types_aggregate_risk.jl 
    CobbDouglasTechnologyAR, HouseholdWorkspace_T12,

# househod_aggregate_risk.jl 
    shadow_rates!,    

# transition_aggregate_risk.jl
    solve_laissez_faire_transition_AR, solve_RPI_transition_AR, solve_RPI_transition_portfolio_AR    

end