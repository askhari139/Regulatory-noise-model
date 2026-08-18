"""
For each of the 10000 RACIPE parameter sets, resample lambda values from the
additive-noise stationary distribution (one draw per sigma level), run 100
deterministic simulations from random ICs, collect steady states.
"""

using DifferentialEquations
using DataFrames
using Statistics
using Random
using ProgressMeter
using CSV

# =============================================================================
# Lambda utilities
# =============================================================================

"""
    get_lambda_indices(param_names)

Return (lambda_indices, max_lambdas, min_lambdas) for all fold-change
parameters found in `param_names`.

Inhibitory (Inh_of_*): lambda in (0.001, 0.999)
Activating (Act_of_*): lambda in (1.001, 100.0)
"""
function get_lambda_indices(param_names::Vector{String})
    lambda_indices = Int[]
    max_lambdas    = Float64[]
    min_lambdas    = Float64[]

    for (i, name) in enumerate(param_names)
        if occursin(r"^Inh_of_", name)
            push!(lambda_indices, i)
            push!(max_lambdas,    0.999)
            push!(min_lambdas,    0.001)
        elseif occursin(r"^Act_of_", name)
            push!(lambda_indices, i)
            push!(max_lambdas,    100.0)
            push!(min_lambdas,    1.001)
        end
    end

    isempty(lambda_indices) && @warn "No fold-change parameters found in param_names"
    return lambda_indices, max_lambdas, min_lambdas
end

"""
    get_threshold_indices(param_names)

Return (threshold_indices, max_thresholds, min_thresholds) for all
threshold parameters (Trd_of_*) in `param_names`.
"""
function get_threshold_indices(param_names::Vector{String})
    trd_indices = Int[]
    max_vals    = Float64[]
    min_vals    = Float64[]

    for (i, name) in enumerate(param_names)
        if occursin(r"^Trd_of_", name)
            push!(trd_indices, i)
            push!(max_vals, 100.0)
            push!(min_vals,   0.01)
        end
    end

    isempty(trd_indices) && @warn "No threshold parameters found in param_names"
    return trd_indices, max_vals, min_vals
end

# =============================================================================
# Core helpers
# =============================================================================

"""
    build_lambda_dist(lambda0, sigma, min_lambda, max_lambda; n_steps, burn_in)

Simulate the additive-noise lambda random walk and return stationary samples
(post burn-in). sigma=0 returns [lambda0].
"""
function build_lambda_dist(lambda0::Float64, sigma::Float64,
                           min_lambda::Float64, max_lambda::Float64;
                           n_steps::Int=50_000, burn_in::Int=10_000)
    if sigma == 0.0
        return rand(n_steps - burn_in) .* (max_lambda - min_lambda) .+ min_lambda
    end

    traj = Vector{Float64}(undef, n_steps)
    lam  = lambda0
    for i in 1:n_steps
        lam     = clamp(lam + randn() * sigma, min_lambda, max_lambda)
        traj[i] = lam
    end
    return traj[burn_in+1:end]
end

"""
    sample_param_vector(base_params, lambda_indices, dists)

Copy of `base_params` with each lambda replaced by one draw from its
pre-computed stationary distribution.
"""
function sample_param_vector(base_params::Vector{Float64},
                             lambda_indices::Vector{Int},
                             dists::Vector{Vector{Float64}})
    p = copy(base_params)
    for (k, idx) in enumerate(lambda_indices)
        p[idx] = rand(dists[k])
    end
    return p
end

"""
    random_ic(param_vector, node_names, param_names)

Random IC scaled to each node's production/degradation ratio.
"""
function random_ic(param_vector::Vector{Float64},
                   node_names::Vector{String},
                   param_names::Vector{String})
    u0 = Vector{Float64}(undef, length(node_names))
    for (j, node) in enumerate(node_names)
        prod_idx = findfirst(x -> x == "Prod_of_$node", param_names)
        deg_idx  = findfirst(x -> x == "Deg_of_$node",  param_names)
        u0[j]    = rand() * 1.5 * param_vector[prod_idx] / param_vector[deg_idx]
    end
    return u0
end

"""
    collect_steady_state(ode_func!, u0, params; tspan, reltol, abstol)

Run a deterministic ODE and return the final state as a proxy for the
steady state.
"""
function collect_steady_state(ode_func!,
                               u0::Vector{Float64},
                               params::Vector{Float64};
                               tspan::Tuple{Float64,Float64}=(0.0, 500.0),
                               reltol::Float64=1e-6,
                               abstol::Float64=1e-7,
                               ss_tol::Float64=1e-2,
                               check_interval::Float64=100.0)
    # Terminate early once all du/dt are below tolerance
    du = similar(u0)
    function steady_state_check!(integrator)
        ode_func!(du, integrator.u, integrator.p, integrator.t)
        if all(abs.(du) .< ss_tol)
            terminate!(integrator)
        end
    end
    cb = PeriodicCallback(steady_state_check!, check_interval)

    prob = ODEProblem(ode_func!, u0, tspan, params)
    sol  = solve(prob, Tsit5(), callback=cb, saveat=tspan[2],
                 reltol=reltol, abstol=abstol)
    return sol[:, end]
end

# =============================================================================
# Per-parameter analysis (called by each SLURM task)
# =============================================================================

"""
    analyze_single_param(ode_func!, param_id, params_data, prs, node_names;
                         sigma_values, num_sims, tspan, seed)

Run the full resampling analysis for one parameter set and return a DataFrame.
"""
function analyze_single_param(ode_func!,
                               param_id::Int,
                               params_data,
                               prs,
                               node_names::Vector{String};
                               sigma_values::Vector{Float64}=[0.0, 0.01, 0.05, 0.1],
                               num_sims::Int=100,
                               tspan::Tuple{Float64,Float64}=(0.0, 500.0),
                               seed::Int=42, 
                               resample::Bool=true)

    Random.seed!(seed + param_id)   # reproducible but distinct per param
    if !resample
        @info "Resampling disabled, using base parameters for ParamID $param_id"
        sigma_values = [0.0]  # only one simulation with base parameters
    end
    param_names    = prs.names
    lambda_indices, max_lambdas, min_lambdas = get_lambda_indices(param_names)

    param_row   = params_data.data[params_data.data.ParamID .== param_id, :]
    nrow(param_row) == 0 && error("ParamID $param_id not found")
    base_params = [param_row[1, Symbol(name)] for name in param_names]

    n_rows    = length(sigma_values) * num_sims
    col_ids   = fill(param_id, n_rows)
    col_sigma = Vector{Float64}(undef, n_rows)
    col_sim   = Vector{Int}(undef, n_rows)
    node_cols = [Vector{Float64}(undef, n_rows) for _ in node_names]

    row = 0
    parList = Vector{Vector{Float64}}() # cache for lambda distributions
    for sigma in sigma_values
        if resample
            dists = [build_lambda_dist(base_params[lambda_indices[k]], sigma,
                                    min_lambdas[k], max_lambdas[k])
                    for k in eachindex(lambda_indices)]
            p  = sample_param_vector(base_params, lambda_indices, dists)
            push!(parList, p)  # cache for lambda distributions
        else
            p = copy(base_params)
            parList = [base_params]  # just one set of parameters, no resampling
        end
        for sim in 1:num_sims
            u0 = random_ic(p, node_names, param_names)
            ss = collect_steady_state(ode_func!, u0, p; tspan=tspan)

            row += 1
            col_sigma[row] = sigma
            col_sim[row]   = sim
            for (j, _) in enumerate(node_names)
                node_cols[j][row] = ss[j]
            end
        end
    end
    parDf = DataFrame(hcat(parList...)', Symbol.(param_names))
    parDf[!, :ParamID] = fill(param_id, nrow(parDf))
    parDf[!, :Sigma]   = sigma_values


    df = DataFrame(ParamID=col_ids, NoiseLevel=col_sigma, SimID=col_sim)
    for (j, name) in enumerate(node_names)
        df[!, Symbol(name)] = node_cols[j]
    end
    return df, parDf
end

# =============================================================================
# SLURM submission
# =============================================================================

"""
    submit_resampled_slurm(network, param_ids, sigma_values, work_dir;
                           num_sims, tspan, log_dir,
                           mem, time_limit, cpus)

Write a param-ID list and a SLURM job-array script, then submit it.
Each array task calls `det_resample_single.jl` for one parameter ID.
"""
function submit_resampled_slurm(network::String,
                                param_ids::Vector{Int},
                                sigma_values::Vector{Float64},
                                work_dir::String;
                                num_sims::Int=100,
                                tspan::Tuple{Float64,Float64}=(0.0, 500.0),
                                log_dir::String=joinpath(work_dir, "logs"),
                                mem::String="4G",
                                time_limit::String="1:00:00",
                                cpus::Int=1, resample::Bool=true)

    mkpath(work_dir)
    mkpath(log_dir)

    sigma_str  = join(sigma_values, ",")
    tspan_str  = "$(tspan[1]),$(tspan[2])"
    ids_file   = joinpath(work_dir, "param_ids.txt")
    source_dir  = get(ENV, "RACIPE_SOURCE", pwd())

    open(ids_file, "w") do f
        for id in param_ids
            println(f, id)
        end
    end

    job_script = """
#!/bin/bash
#SBATCH --job-name=det_resample_$(network)
#SBATCH --output=$(log_dir)/param_%a.out
#SBATCH --error=$(log_dir)/param_%a.err
#SBATCH --array=1-$(length(param_ids))
#SBATCH --cpus-per-task=$(cpus)
#SBATCH --mem=$(mem)
#SBATCH --time=$(time_limit)
#SBATCH --partition=ctbp

cd $(pwd())

PARAM_ID=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $(ids_file))

echo "Task \${SLURM_ARRAY_TASK_ID}: ParamID \${PARAM_ID}"

julia --project=. --threads=$(cpus) scripts/det_resample_single.jl network=$(network) param-id=\${PARAM_ID} sigma-values=$(sigma_str) work-dir=$(work_dir) num-sims=$(num_sims) tspan=$(tspan_str) resample=$(resample)
"""

    # Detect cluster max array size
    max_array_size = 1000  # conservative default
    try
        config_out = read(`scontrol show config`, String)
        m = match(r"MaxArraySize\s*=\s*(\d+)", config_out)
        if m !== nothing
            max_array_size = parse(Int, m.captures[1]) - 1  # SLURM upper bound is exclusive
        end
    catch
        @warn "Could not query MaxArraySize, defaulting to $max_array_size"
    end

    n_total  = length(param_ids)
    n_chunks = ceil(Int, n_total / max_array_size)

    println("Submitting $n_total tasks in $n_chunks array(s) (max size: $max_array_size)")

    for chunk in 1:n_chunks
        i_start  = (chunk - 1) * max_array_size + 1
        i_end    = min(chunk * max_array_size, n_total)
        chunk_size = i_end - i_start + 1

        chunk_script = """
#!/bin/bash
#SBATCH --job-name=det_resample_$(network)_$(chunk)
#SBATCH --output=$(log_dir)/param_%a.out
#SBATCH --error=$(log_dir)/param_%a.err
#SBATCH --array=1-$(chunk_size)
#SBATCH --cpus-per-task=$(cpus)
#SBATCH --mem=$(mem)
#SBATCH --time=$(time_limit)
#SBATCH --partition=ctbp

cd $(pwd())

# Offset array task ID into the global param_ids list
GLOBAL_IDX=\$(( \${SLURM_ARRAY_TASK_ID} + $(i_start - 1) ))
PARAM_ID=\$(sed -n "\${GLOBAL_IDX}p" $(ids_file))

echo "Task \${SLURM_ARRAY_TASK_ID} (global \${GLOBAL_IDX}): ParamID \${PARAM_ID}"

julia --project=. --threads=$(cpus) scripts/det_resample_single.jl network=$(network) param-id=\${PARAM_ID} sigma-values=$(sigma_str) work-dir=$(work_dir) num-sims=$(num_sims) tspan=$(tspan_str) resample=$resample
"""
        chunk_file = joinpath(work_dir, "submit_array_$(chunk).sh")
        open(chunk_file, "w") do f
            write(f, chunk_script)
        end

        run(`sbatch $(chunk_file)`)
        println("  Submitted chunk $chunk/$n_chunks (tasks $i_start – $i_end)")
    end

    println("Done. Monitor with: squeue -u \$USER")
    println("Logs: $(log_dir)/")
end