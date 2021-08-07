module Aiyagari

using UnPack
using StaticArrays
using StructArrays
using Roots
using ProgressMeter
using LoopVectorization
using LinearAlgebra
using FastClosures
using FLoops
using Tullio
using NLsolve
using Plots
using StatsBase
using LaTeXStrings


include("types.jl")
include("calibration.jl")
include("helpers.jl")
include("firms.jl")
include("fiscal.jl")
include("preferences.jl")
include("household_problem.jl")
include("distribution_and_aggregates.jl")
include("stationary_equilibrium.jl")
include("transition.jl")
include("summaries_plots.jl")


const _TOL_VALUE = 1e-6
const _TOL_PDF = 1e-8


export

# types.jl

    CobbDouglas,
    CES,
    Technology,
    MarkupTechnology,

    LinearIncomeTaxes,

    CRRA,
    EZ,
    GHH,
    FixedLabor,

    Household,

    Economy,

    get_e,
    get_h,
    get_t,

# household_problem.jl

    solve_stationary_household!,
    solve_stationary_household,
    is_valid,
    consumption_alloc,

# distribution_and_aggregates.jl

    stationary_pdf!,
    stationary_pdf,
    asset_supply,
    labor_supply,
    aggregate,

# firms.jl

    output,
    mpk_from_factors,
    golden_rule_k,

    rK_from_r,
    mpk_from_after_tax_rK,
    k_from_mpk,

# fiscal.jl

    taxes,

# calibration.jl

    calibration,

# stationary_equilibrium.jl

    solve_laissez_faire,
    solve_new_stationary_equilibrium_given_k_b,

# transition.jl

    solve_transition,

# summaries_plots.jl

    summary_statics,
    do_plots

end
