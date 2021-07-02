module Aiyagari 

using UnPack
using StaticArrays
using StructArrays
using Roots
using ProgressMeter
using LoopVectorization
using LinearAlgebra
using FastClosures
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
include("household.jl") 
include("stationaryequilibrium.jl")
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

# households.jl

    solve_stationary_household!,
    solve_stationary_household,

    stationary_pdf!, 
    stationary_pdf, 

    asset_supply, 
    labor_supply,

    is_valid,

    consumption,
    aggregate,

# firms.jl

    get_y,
    get_mpk,
    get_golden_k,

    rK_from_r,
    mpk_from_after_tax_rK,
    k_from_mpk,

# fiscal.jl

    get_taxes, 

# calibration.jl

    calibration,

# stationaryequilibrium.jl

    solve_laissez_faire,
    solve_new_stationary_equilibrium_given_k_b, 

# transition.jl

    solve_transition,

# summaries_plots.jl

    summary_statics,
    do_plots
    
end 
