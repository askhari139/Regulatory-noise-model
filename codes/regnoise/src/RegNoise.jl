"""
RegNoise — stochastic regulatory-strength (λ) simulations on RACIPE networks.
See SPEC.md for the full specification; every number comes from the config file.
"""
module RegNoise

using CSV, DataFrames, JLD2, TOML, SHA, Dates, Random, Statistics, StatsBase
using OrdinaryDiffEq, DiffEqCallbacks, SciMLBase

include("config.jl")
include("racipe.jl")
include("simulate.jl")
include("observables.jl")
include("conditions.jl")
include("io.jl")
include("selection.jl")

export load_config, verify_inputs, load_network, select_parameters, load_selection,
       read_parameters, run_unit, write_unit, unit_status, validate_unit, raw_path,
       lambda_bounds, check_clean_tree, run_dir, expected_cells

end
