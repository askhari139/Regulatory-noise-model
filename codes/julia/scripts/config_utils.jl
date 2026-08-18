"""
config_utils.jl

Load/save helpers for pipeline configuration files.
Uses JLD2 so all Julia types (Dict{Int,String}, Symbol, Vector{Float64}, …)
round-trip without any manual conversion.
"""

using JLD2
using JSON
# ============================================================
# PARAMETER CONFIG  (written by 1_select_parameters.jl)
# ============================================================

"""
    save_parameter_config(path, data::Dict)

Save parameter-selection results to a JLD2 file.

Expected keys in `data`:
  "network"           => String
  "param_ids"         => Vector{Int}
  "param_types"       => Dict{Int,String}
  "node_names"        => Vector{String}
  "lambda_indices"    => Vector{Int}
  "max_lambdas"       => Vector{Float64}
  "min_lambdas"       => Vector{Float64}
  "racipe_thresholds" => Vector{Float64}
  "param_type"        => String  (or Symbol — both accepted)
  "timestamp"         => String  (optional)
"""
function save_parameter_config(path::String, data::Dict)
    jldsave(path;
        network           = String(data["network"]),
        param_ids         = Vector{Int}(data["param_ids"]),
        param_types       = Dict{Int,String}(data["param_types"]),
        node_names        = Vector{String}(data["node_names"]),
        lambda_indices    = Vector{Int}(data["lambda_indices"]),
        max_lambdas       = Vector{Float64}(data["max_lambdas"]),
        min_lambdas       = Vector{Float64}(data["min_lambdas"]),
        racipe_thresholds = Vector{Float64}(data["racipe_thresholds"]),
        param_type        = String(data["param_type"]),
        timestamp         = string(get(data, "timestamp", "")),
        data_folder       = data["data_folder"],
    )
end

"""
    load_parameter_config(path) -> Dict

Load a parameter config saved by `save_parameter_config`.
All values are returned with their native Julia types — no string-to-Int
key conversion needed.
"""
function load_parameter_config(path::String)
    d = load(path)
    return Dict{String,Any}(
        "network"           => d["network"]::String,
        "param_ids"         => d["param_ids"]::Vector{Int},
        "param_types"       => d["param_types"]::Dict{Int,String},
        "node_names"        => d["node_names"]::Vector{String},
        "lambda_indices"    => d["lambda_indices"]::Vector{Int},
        "max_lambdas"       => d["max_lambdas"]::Vector{Float64},
        "min_lambdas"       => d["min_lambdas"]::Vector{Float64},
        "racipe_thresholds" => d["racipe_thresholds"]::Vector{Float64},
        "param_type"        => d["param_type"]::String,
        "timestamp"         => get(d, "timestamp", ""),
        "data_folder"       => d["data_folder"],
    )
end

# ============================================================
# JOB CONFIG  (written by 2_submit_jobs.jl)
# ============================================================

"""
    save_job_config(path, data::Dict)

Save job submission configuration to a JLD2 file.

Expected keys in `data`:
  "network"      => String
  "noise_mode"   => String
  "param_type"   => String
  "num_sims"     => Int
  "noise_levels" => Vector{Float64}
  "param_ids"    => Vector{Int}
  "param_types"  => Dict{Int,String}
  "timestamp"    => String  (optional)
"""
function save_job_config(path::String, data::Dict)
    jldsave(path;
        network      = String(data["network"]),
        noise_mode   = String(data["noise_mode"]),
        param_type   = String(data["param_type"]),
        num_sims     = Int(data["num_sims"]),
        noise_levels = Vector{Float64}(data["noise_levels"]),
        param_ids    = Vector{Int}(data["param_ids"]),
        param_types  = Dict{Int,String}(data["param_types"]),
        timestamp    = string(get(data, "timestamp", "")),
        dt           = data["dt"],

    )
end

"""
    load_job_config(path) -> Dict

Load a job config saved by `save_job_config`.
"""
function load_job_config(path::String)
    d = load(path)
    return Dict{String,Any}(
        "network"      => d["network"]::String,
        "noise_mode"   => d["noise_mode"]::String,
        "param_type"   => d["param_type"]::String,
        "num_sims"     => d["num_sims"]::Int,
        "noise_levels" => d["noise_levels"]::Vector{Float64},
        "param_ids"    => d["param_ids"]::Vector{Int},
        "param_types"  => d["param_types"]::Dict{Int,String},
        "timestamp"    => get(d, "timestamp", ""),
        "dt"           => d["dt"],
    )
end

# ============================================================
# MIGRATION HELPER  (one-time use)
# ============================================================

"""
    migrate_json_to_jld2(json_path; delete_json=false)

Convert an existing JSON config file produced by the old pipeline to JLD2.
Pass the path to the `.json` file; the `.jld2` file is written alongside it.

Set `delete_json=true` to remove the original after a successful conversion.

Example:
    migrate_json_to_jld2("Additive/TS/selected_parameters.json")
    migrate_json_to_jld2("Additive/TS/job_config.json")
"""
function migrate_json_to_jld2(json_path::String; delete_json::Bool=false)
      # lazy import — only needed during migration

    jld2_path = replace(json_path, r"\.json$" => ".jld2")

    raw = open(JSON.parse, json_path)

    if occursin("selected_parameters", json_path)
        # Reconstruct with correct types
        data = Dict{String,Any}(
            "network"           => String(raw["network"]),
            "param_ids"         => Vector{Int}(raw["param_ids"]),
            "param_types"       => Dict{Int,String}(
                                       parse(Int, k) => String(v)
                                       for (k, v) in raw["param_types"]),
            "node_names"        => Vector{String}(raw["node_names"]),
            "lambda_indices"    => Vector{Int}(raw["lambda_indices"]),
            "max_lambdas"       => Vector{Float64}(raw["max_lambdas"]),
            "min_lambdas"       => Vector{Float64}(raw["min_lambdas"]),
            "racipe_thresholds" => Vector{Float64}(raw["racipe_thresholds"]),
            "param_type"        => String(raw["param_type"]),
            "timestamp"         => string(get(raw, "timestamp", "")),
        )
        save_parameter_config(jld2_path, data)

    elseif occursin("job_config", json_path)
        data = Dict{String,Any}(
            "network"      => String(raw["network"]),
            "noise_mode"   => String(raw["noise_mode"]),
            "param_type"   => String(raw["param_type"]),
            "num_sims"     => Int(raw["num_sims"]),
            "noise_levels" => Vector{Float64}(raw["noise_levels"]),
            "param_ids"    => Vector{Int}(raw["param_ids"]),
            "param_types"  => Dict{Int,String}(
                                  parse(Int, k) => String(v)
                                  for (k, v) in raw["param_types"]),
            "timestamp"    => string(get(raw, "timestamp", "")),
        )
        save_job_config(jld2_path, data)

    else
        error("Unrecognised config file: $json_path\n" *
              "Expected a path containing 'selected_parameters' or 'job_config'.")
    end

    println("✓ Migrated  $json_path  →  $jld2_path")
    if delete_json
        rm(json_path)
        println("  (deleted original)")
    end
    return jld2_path
end