"""
Stochastic simulations with parameter noise for IDP modeling
Implements the methodology from the paper for toggle switch and general networks
"""
module StochasticSimulations

using OrdinaryDiffEq
using DiffEqCallbacks
using Distributions
using Statistics
using DataFrames
using Plots
using Random
using Distributed
import StatsBase.skewness, StatsBase.kurtosis

export simulate_with_noise, simulate_with_noise_tracked, run_multiple_stochastic_simulations, create_summary_dataframe
export discretize_trajectory, calculate_mrt, count_switches, aggregate_mrt, aggregate_transitions
export StochasticResult, analyze_noise_effects, analyze_noise_effects_parallel, analyze_noise_effects_distributed, analyze_noise_effects_det
# export plot_trajectories, plot_mrt_comparison, plot_switches_vs_noise
export get_lambda_indices_from_topology, get_lambda_dist
# export plot_lambda_trajectories, plot_lambda_distributions
# export plot_lambda_vs_states, plot_lambda_summary
# ============================================================
# Data Structures
# ============================================================
"""
    StochasticResult

Container for results of stochastic simulation
"""
struct StochasticResult
    # sol::ODESolution                    # Full ODE solution
    # times::Vector{Float64}              # Time points
    # states::Matrix{Float64}             # State trajectories [nodes × time]
    discrete_states::Vector{<:Tuple}      # Discretized states
    thresholds::Vector{Float64}         # Thresholds for discretization
    mrt::Dict{Tuple, Float64}          # Mean residence time
    switches::Int                       # Number of switching events
    transition_counts::Dict
    noise_level::Float64               # sigma used
    trajectory_stats::Dict{String, Vector{Float64}}   # stat name => per-node values
end

# ============================================================
# Core Simulation Functions
# ============================================================

"""
    simulate_with_noise(ode_func!, u0, tspan, params, lambda_indices; 
                       sigma=0.01, dt=0.01, saveat=0.01, 
                       keep_positive=true, min_lambda=0.001)

Run ODE simulation with Gaussian noise on specified parameters

# Arguments
- `ode_func!`: ODE function (du, u, p, t) -> nothing
- `u0`: Initial conditions
- `tspan`: Time span (t_start, t_end)
- `params`: Parameter vector
- `lambda_indices`: Indices of lambda parameters to perturb
- `sigma`: Standard deviation of noise (default: 0.01)
- `dt`: Time step for noise application (default: 0.01)
- `saveat`: Save interval for solution (default: 0.01)
- `keep_positive`: Keep lambda > min_lambda (default: true)
- `min_lambda`: Minimum value for lambda (default: 0.001)

# Returns
- ODESolution object with stochastic parameter trajectories
"""
function simulate_with_noise(ode_func!, u0, tspan, params, lambda_indices;
                            sigma::Float64=0.01,
                            dt::Float64=0.01,
                            saveat::Float64=0.01,
                            keep_positive::Bool=true,
                            min_lambdas::Vector{Float64}=Float64[],
                            max_lambdas::Vector{Float64}=Float64[],
                            noise_mode::String="Additive")
    
    p_current = copy(params)
    original_params = params[lambda_indices]
 
    if sigma == 0.0
        prob = ODEProblem(ode_func!, u0, tspan, p_current)
        return solve(prob, Tsit5(), saveat=saveat, reltol=1e-4, abstol=1e-5)
    end
 
    # Fix: don't shadow `l`, and build clean local copies
    local_max_lambdas = copy(max_lambdas)
    local_min_lambdas = copy(min_lambdas)
 
    if length(local_max_lambdas) < length(lambda_indices)
        n_missing = length(lambda_indices) - length(local_max_lambdas)
        append!(local_max_lambdas, fill(Inf, n_missing))
        for (i, idx) in enumerate(lambda_indices)       # ← `idx` not `l`
            if p_current[idx] < 1
                local_max_lambdas[i] = 0.999
            end
        end
    end
 
    if length(local_min_lambdas) < length(lambda_indices)
        n_missing = length(lambda_indices) - length(local_min_lambdas)
        append!(local_min_lambdas, fill(0.001, n_missing))
        for (i, idx) in enumerate(lambda_indices)       # ← `idx` not `l`
            if p_current[idx] > 1
                local_min_lambdas[i] = 1.001
            end
        end
    end
 
    function affect!(integrator)
        for (i, idx) in enumerate(lambda_indices)
            effective_sigma = local_max_lambdas[i] > 1.0 ? sigma * local_max_lambdas[i] : sigma
            noise = randn() * effective_sigma
            if noise_mode == "Additive"
                integrator.p[idx] += noise
            elseif noise_mode == "Multiplicative"
                integrator.p[idx] *= (1 + noise)
            elseif noise_mode == "MultiplicativeInvLambda"
                # Plain "Multiplicative" scales sigma by the link's static
                # max_lambda bound (100 for activation), which is large relative
                # to lambda values near the weak/no-effect boundary (~1). Since
                # lambda*(1+noise) then clamps at that nearby lower bound far
                # more often than at the distant upper bound, activation links
                # get systematically eroded back toward "no effect" instead of
                # perturbed symmetrically. Scaling sigma by max_lambda/lambda_now
                # (inversely with the CURRENT lambda, not the static bound) keeps
                # the kick large near baseline (same as today) but shrinks it
                # once a link has already strengthened, so strong activators
                # stop being an easy target for the next kick back down.
                # Inhibition links (max_lambda <= 1) are unaffected.
                eff_sigma_inv = local_max_lambdas[i] > 1.0 ?
                    sigma * local_max_lambdas[i] / integrator.p[idx] : sigma
                integrator.p[idx] *= (1 + randn() * eff_sigma_inv)
            elseif noise_mode == "Lognormal"
                integrator.p[idx] *= exp(noise)
            elseif noise_mode == "Fluctuating"
                integrator.p[idx] = original_params[i] + noise
            elseif noise_mode == "Jumping"
                integrator.p[idx] += (noise < 0 ? -effective_sigma : effective_sigma)
            elseif noise_mode == "Extreme"
                integrator.p[idx] += (noise < 0 ? -1.0 : 1.0)
            else
                integrator.p[idx] += noise
            end
 
            # Clamp to bounds
            integrator.p[idx] = clamp(integrator.p[idx], local_min_lambdas[i], local_max_lambdas[i])
        end
    end
    
    cb = PeriodicCallback(affect!, dt; save_positions = (false, false))
    prob = ODEProblem(ode_func!, u0, tspan, p_current)
    
    if noise_mode in ("Jumping", "Extreme")
        return solve(prob, Tsit5(), callback=cb, saveat=saveat,
                     tstops=collect(tspan[1]:dt:tspan[2]),
                     reltol=1e-4, abstol=1e-5)
    else
        return solve(prob, Tsit5(), callback=cb, saveat=saveat,
                     reltol=1e-4, abstol=1e-5)
    end
end


function simulate_steady_state(ode_func!, u0, tspan, params;
                               saveat::Float64=0.01,
                               n_tail::Int=100)

    prob = ODEProblem(ode_func!, u0, tspan, params)
    sol  = solve(prob, Tsit5(), saveat=saveat, reltol=1e-4, abstol=1e-5)

    n_pts   = length(sol.t)
    i_start = max(1, n_pts - n_tail + 1)

    # Return a plain matrix (nodes × n_tail) and the corresponding time points
    return Array(sol)[:, i_start:end], sol.t[i_start:end]
end


"""
    simulate_with_noise_tracked(ode_func!, u0, tspan, params, lambda_indices; 
                                sigma=0.01, dt=0.01, saveat=0.01)

Simulate with noise and track lambda parameter values over time

# Returns
- `sol`: ODE solution
- `lambda_history`: Matrix of lambda values over time (length(lambda_indices) × n_timepoints)
- `times`: Time points where lambda was recorded
"""
function simulate_with_noise_tracked(ode_func!, u0, tspan, params, lambda_indices;
                                    sigma::Float64=0.01,
                                    dt::Float64=0.01,
                                    saveat::Float64=0.01,
                                    keep_positive::Bool=true,
                                    min_lambdas::Vector{Float64}=Float64[],
                                    max_lambdas::Vector{Float64}=Float64[],
                            noise_mode::String="Additive")
    
    # Create mutable copy of parameters
    p_current = copy(params)
    original_params = params[lambda_indices]

    # Storage for lambda history
    lambda_history = Vector{Vector{Float64}}()
    time_history = Float64[]
    if length(max_lambdas) < length(lambda_indices)
        l = length(lambda_indices) - length(max_lambdas)
        max_lambdas = vcat(max_lambdas, fill(Inf, l))
        for (i,l) in enumerate(lambda_indices)
            if p_current[l] < 1
                max_lambdas[i] = 0.999
            end
        end
    end
    if length(min_lambdas) < length(lambda_indices)
        l = length(lambda_indices) - length(min_lambdas)
        min_lambdas = vcat(min_lambdas, fill(0.001, l))
        for (i,l) in enumerate(lambda_indices)
            if p_current[l] > 1
                min_lambdas[i] = 1.001
            end
        end
    end
    # Callback to add noise and record lambda values
    function affect!(integrator)
        # Record current time and lambda values
        push!(time_history, integrator.t)
        push!(lambda_history, [integrator.p[idx] for idx in lambda_indices])
        
        # Add noise
        for (i,idx) in enumerate(lambda_indices)
            noise = rand(Normal(0.0, sigma))
            if noise_mode == "Additive"
                integrator.p[idx] += noise
            elseif noise_mode == "Multiplicative"
                integrator.p[idx] *= (1 + noise)
            elseif noise_mode == "Lognormal"
                integrator.p[idx] *= exp(noise)
            elseif noise_mode == "Fluctuating"
                integrator.p[idx] = original_params[i] + noise
            else
                integrator.p[idx] += noise
            end

            
            # Keep parameters positive if requested
            if keep_positive && integrator.p[idx] < min_lambdas[i]
                integrator.p[idx] = min_lambdas[i]
            end
            if integrator.p[idx] > max_lambdas[i]
                integrator.p[idx] = max_lambdas[i]
            end
        end
    end
    
    # Setup periodic callback
    cb = PeriodicCallback(affect!, dt)
    
    # Solve with callback
    prob = ODEProblem(ode_func!, u0, tspan, p_current)
    sol = solve(prob, Tsit5(), callback=cb, saveat=saveat)
    
    # Convert lambda_history to matrix (n_lambda × n_timepoints)
    lambda_matrix = hcat(lambda_history...)
    
    return sol, lambda_matrix, time_history
end

"""
    get_lambda_indices_from_topology(param_names, topo)

Automatically identify lambda (fold-change) parameter indices from parameter names
For inhibitions: Inh_of_X, for activations: Act_of_X

# Arguments
- `param_names`: Vector of parameter names (from RACIPE .prs file)
- `topo`: Topology object (optional, for validation)

# Returns
- Vector of indices corresponding to fold-change parameters
"""
function get_lambda_indices_from_topology(param_names::Vector{String}; param_type::Symbol = :FoldChange)
    lambda_indices = Int[]
    max_lambdas = Float64[]
    min_lambdas = Float64[]
    for (i, name) in enumerate(param_names)
        # Look for inhibition or activation fold-change parameters
        if param_type == :FoldChange
            if occursin(r"^Inh_of_", name)
                push!(lambda_indices, i)
                push!(max_lambdas, 0.999)
                push!(min_lambdas, 0.001)
            end
            if occursin(r"^Act_of_", name)
                push!(lambda_indices, i)
                push!(max_lambdas, 100.0)
                push!(min_lambdas, 1.001)
            end
        elseif param_type == :Threshold
            if occursin(r"^Trd_of_", name)
                push!(lambda_indices, i)
                push!(max_lambdas, 100)
                push!(min_lambdas, 0.01)
            end
        end
    end
    
    if isempty(lambda_indices)
        @warn "No fold-change (lambda) parameters found in parameter names"
    end
    
    return lambda_indices, max_lambdas, min_lambdas
end

function get_lambda_dist(noise_mode::String, noise_level::Float64, min_lambda::Float64, max_lambda::Float64)
    lambda_init = collect(1:100)*1.0
    if min_lambda < 1
        lambda_init = lambda_init./100
    end
    # lambda_init_key = round.(lambda_init, digits = 0)
    lambda_trajectories = Dict{Float64, Vector{Float64}}()
    for i in eachindex(lambda_init)
        lambda = lambda_init[i]
        key = copy(lambda_init[i])
        lt = zeros(10000)
        trajectories = Float64[]
        for iter in 1:100
            noise = noise_level*randn(10000)
            lt[1] = lambda
            for t in 2:10000
                if noise_mode == "Additive"
                    ln = lt[t-1] + noise[t]
                elseif noise_mode == "Multiplicative"
                    ln = lt[t-1]*(1 + noise[t])
                elseif noise_mode == "Fluctuating"
                    ln = lambda + noise[t]
                elseif noise_mode == "Jumping"
                    ln = lt[t-1] + noise[t] > 0 ? noise_level : -1*noise_level
                elseif noise_mode == "Extreme"
                    ln = lt[t-1] + noise[t] > 0 ? 1 : -1
                end
                lt[t] = max(min(ln, max_lambda), min_lambda)
            end
            trajectories = vcat(trajectories, lt[5000:end])
        end
        lambda_trajectories[key] = trajectories
    end
    return(lambda_trajectories)
end

# ============================================================
# Discretization and Analysis
# ============================================================

"""
    discretize_trajectory(sol, thresholds)

Convert continuous trajectory to discrete states

# Arguments
- `sol`: ODE solution
- `thresholds`: Vector of thresholds for each node

# Returns
- Vector of tuples representing discrete states at each time point
"""
function discretize_trajectory(sol::ODESolution, thresholds::Vector{Float64})
    sol_matrix = Array(sol)  # [nodes × timepoints]
    
    # Broadcast comparison with reshaped thresholds
    threshold_matrix = repeat(thresholds, 1, size(sol_matrix, 2))
    discrete_matrix = Int.(sol_matrix .> threshold_matrix)
    
    # Convert to vector of tuples
    discrete_states = [Tuple(discrete_matrix[:, i]) for i in 1:size(discrete_matrix, 2)]
    
    return discrete_states
end

"""
    calculate_mrt(discrete_states)

Calculate Mean Residence Time for each state

# Arguments
- `discrete_states`: Vector of discrete state tuples

# Returns
- Dictionary mapping state => residence time (fraction of total time)
"""
function calculate_mrt(discrete_states::Vector{<:Tuple})
    total_time = length(discrete_states)
    state_counts = Dict{Tuple, Int}()
    
    # Count occurrences of each state
    for state in discrete_states
        state_counts[state] = get(state_counts, state, 0) + 1
    end
    
    # Convert to fractions
    mrt = Dict{Tuple, Float64}()
    for (state, count) in state_counts
        mrt[state] = count / total_time
    end
    
    return mrt
end

"""
    count_switches(discrete_states)

Count number of state transitions

# Arguments
- `discrete_states`: Vector of discrete state tuples

# Returns
- Number of switching events
"""
function count_switches(discrete_states::Vector{<:Tuple})
    return sum(discrete_states[1:end-1] .!= discrete_states[2:end])
end

function count_transitions(discrete_states::Vector{<:Tuple})
    """
    Count transition events between discrete states
    
    Returns:
    - Dict of (from_state, to_state) => count
    """
    transition_counts = Dict{Tuple{eltype(discrete_states), eltype(discrete_states)}, Int}()
    
    for i in 2:length(discrete_states)
        if discrete_states[i] != discrete_states[i-1]
            transition = (discrete_states[i-1], discrete_states[i])
            transition_counts[transition] = get(transition_counts, transition, 0) + 1
        end
    end
    
    return transition_counts
end

"""
    average_dwell_time(discrete_states)

Calculate average time spent in each state before switching

# Returns
- Dict of state => average dwell time
"""
function average_dwell_time(discrete_states::Vector{<:Tuple})
    dwell_times = Dict{Tuple, Vector{Int}}()
    
    current_state = discrete_states[1]
    current_dwell = 1
    
    for i in 2:length(discrete_states)
        if discrete_states[i] == current_state
            current_dwell += 1
        else
            # State changed
            if !haskey(dwell_times, current_state)
                dwell_times[current_state] = Int[]
            end
            push!(dwell_times[current_state], current_dwell)
            
            current_state = discrete_states[i]
            current_dwell = 1
        end
    end
    
    # Don't forget the last state
    if !haskey(dwell_times, current_state)
        dwell_times[current_state] = Int[]
    end
    push!(dwell_times[current_state], current_dwell)
    
    # Calculate averages
    avg_dwell = Dict{Tuple, Float64}()
    for (state, times) in dwell_times
        avg_dwell[state] = mean(times)
    end
    
    return avg_dwell
end


function get_trajectory_stats(sol::ODESolution; cut_fraction::Float64=0.5)
    sol_matrix = Array(sol)                       # nodes × time
    n_t        = size(sol_matrix, 2)
    start      = max(1, Int(round(cut_fraction * n_t)))
    tail       = sol_matrix[:, start:end]         # nodes × tail_length

    return Dict(
        "Mean"     => [mean(row)     for row in eachrow(tail)],
        "StDev"    => [std(row)      for row in eachrow(tail)],
        "Skewness" => [skewness(row) for row in eachrow(tail)],
        "Kurtosis" => [kurtosis(row) for row in eachrow(tail)],
    )
end
# ============================================================
# High-Level Analysis Functions
# ============================================================

"""
    run_multiple_stochastic_simulations(ode_func!, params, lambda_indices, node_names;
                                       num_sims=25, tspan=(0.0, 1000.0), sigma=0.01,
                                       threshold_method=:from_sims)

Run multiple stochastic simulations and aggregate results

# Arguments
- `ode_func!`: ODE function
- `params`: Parameter vector
- `lambda_indices`: Indices of lambda parameters
- `node_names`: Names of nodes (for output)
- `num_sims`: Number of simulations (default: 25)
- `tspan`: Time span (default: (0.0, 1000.0))
- `sigma`: Noise level (default: 0.01)
- `threshold_method`: :from_sims (calculate from all sims) or :per_sim

# Returns
- Vector of StochasticResult objects
"""
function run_multiple_stochastic_simulations(ode_func!, params, param_names, lambda_indices, node_names;
                                            num_sims::Int=25,
                                            tspan::Tuple{Float64,Float64}=(0.0, 1000.0),
                                            sigma::Float64=0.01,
                                            dt::Float64=0.01,
                                            saveat::Float64=1.0,
                                            u0List::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
                                            racipe_thresholds::Vector{Float64}=Vector{Float64}(),
                                            # threshold_method::Symbol=:from_sims,
                                            max_lambdas::Vector{Float64}=Float64[],
                                            min_lambdas::Vector{Float64}=Float64[],
                                            track_lambda::Bool=false,
                                            seed::Int=123, noise_mode::String = "Additive", 
                                            getSolutions::Bool = false,
                                            cut_fraction::Float64 = 0.5)
    
    # Random.seed!(seed)
    
    n_nodes = length(node_names)
    results = StochasticResult[]
    
    # First, run all simulations
    # println("Running $num_sims simulations with sigma=$sigma...")
    solutions = []
    if isempty(u0List)
        randInit = true
    else
        compatibles = findall(x -> length(x) == length(node_names), u0List)
        if isempty(compatibles)
            randInit = true
        else
            u0List = u0List[compatibles]
            randInit = false
        end
    end

    if length(u0List) < num_sims
        u0List = u0List[mod1.(1:num_sims, length(u0List))]
    end
    
    for i in 1:num_sims
        # Random initial condition
        # Between 0 and 1.5 * steady state approximation (as in paper)
        if randInit
            u0 = zeros(n_nodes)
            for j in 1:n_nodes
                # Rough approximation: steady state ~ production/degradation
                prod_idx = findfirst(x -> occursin("Prod_of_$(node_names[j])", x), 
                                    param_names)
                deg_idx = findfirst(x -> occursin("Deg_of_$(node_names[j])", x),
                                param_names)
                
                # Simple random initial condition
                u0[j] = rand() * 1.5* params[prod_idx]/params[deg_idx]
            end
        else
            u0 = u0List[i]
        end
        
        # Simulate
        if track_lambda
            sol = simulate_with_noise_tracked(ode_func!, u0, tspan, params, lambda_indices,
                                sigma=sigma, dt=dt, saveat=saveat, max_lambdas=max_lambdas, min_lambdas=min_lambdas,noise_mode=noise_mode)
        else
            sol = simulate_with_noise(ode_func!, u0, tspan, params, lambda_indices,
                                sigma=sigma, dt=dt, saveat=saveat, max_lambdas=max_lambdas, min_lambdas=min_lambdas,noise_mode = noise_mode)
        end
        push!(solutions, sol)
        
        # if i % 5 == 0
        #     print(".")
        # end
    end
    # println(" Done!")
    
    # Calculate thresholds
    if length(racipe_thresholds) == length(node_names)
        thresholds = copy(racipe_thresholds)
    else
        @warn "Improper threshodls. Cannot discretize. Exiting..."
        return
    end
    
    # Analyze each simulation
    # println("\nAnalyzing results...")
    for (i, sol) in enumerate(solutions)
        # Discretize
        discrete_states = discretize_trajectory(sol, thresholds)
        n_cut = floor(Int, (1 - cut_fraction) * length(discrete_states))
        discrete_states = discrete_states[n_cut+1:end]
        stats = get_trajectory_stats(sol; cut_fraction = cut_fraction)
        # Calculate metrics
        mrt = calculate_mrt(discrete_states)
        switches = count_switches(discrete_states)
        transitions = count_transitions(discrete_states)
        
        # Store result
        result = StochasticResult(
            # sol,
            # collect(sol.t),
            # Array(sol),
            discrete_states,
            thresholds,
            mrt,
            switches,
            transitions,
            sigma,
            stats
        )
        
        push!(results, result)
    end
    if !getSolutions
        return results
    else
        # Pre-allocate by building all rows at once
        all_times = Float64[]
        all_states = [Float64[] for _ in node_names]
        all_sim_ids = Int[]

        for (i, sol) in enumerate(solutions)
            t = collect(sol.t)
            states = Array(sol)  # nodes × timepoints
            append!(all_times, t)
            append!(all_sim_ids, fill(i, length(t)))
            for (j, _) in enumerate(node_names)
                append!(all_states[j], states[j, :])
            end
        end

        sol_df = DataFrame(SimID = all_sim_ids, Time = all_times)
        for (j, node) in enumerate(node_names)
            sol_df[!, Symbol(node)] = all_states[j]
        end

        return results, sol_df
    end
end


function run_multiple_deterministic_simulations(ode_func!, params, param_names, lambda_indices, node_names;
                                            num_sims::Int=25,
                                            tspan::Tuple{Float64,Float64}=(0.0, 1000.0),
                                            sigma::Float64=0.01,
                                            dt::Float64=0.01,
                                            saveat::Float64=1.0,
                                            u0List::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
                                            racipe_thresholds::Vector{Float64}=Vector{Float64}(),
                                            # threshold_method::Symbol=:from_sims,
                                            max_lambdas::Vector{Float64}=Float64[],
                                            min_lambdas::Vector{Float64}=Float64[],
                                            track_lambda::Bool=false,
                                            seed::Int=123, noise_mode::String = "Additive", 
                                            getSolutions::Bool = false,
                                            cut_fraction::Float64 = 0.5)
    
    # Random.seed!(seed)
    
    n_nodes = length(node_names)
    results = StochasticResult[]
    
    # First, run all simulations
    # println("Running $num_sims simulations with sigma=$sigma...")
    solutions = []
    if isempty(u0List)
        randInit = true
    else
        compatibles = findall(x -> length(x) == length(node_names), u0List)
        if isempty(compatibles)
            randInit = true
        else
            u0List = u0List[compatibles]
            randInit = false
        end
    end

    if length(u0List) < num_sims
        u0List = u0List[mod1.(1:num_sims, length(u0List))]
    end
    
    for i in 1:num_sims
        # Random initial condition
        # Between 0 and 1.5 * steady state approximation (as in paper)
        if randInit
            u0 = zeros(n_nodes)
            for j in 1:n_nodes
                # Rough approximation: steady state ~ production/degradation
                prod_idx = findfirst(x -> occursin("Prod_of_$(node_names[j])", x), 
                                    param_names)
                deg_idx = findfirst(x -> occursin("Deg_of_$(node_names[j])", x),
                                param_names)
                
                # Simple random initial condition
                u0[j] = rand() * 1.5* params[prod_idx]/params[deg_idx]
            end
        else
            u0 = u0List[i]
        end
        
        # Simulate
        # if track_lambda
        #     sol = simulate_with_noise_tracked(ode_func!, u0, tspan, params, lambda_indices,
        #                         sigma=sigma, dt=dt, saveat=saveat, max_lambdas=max_lambdas, min_lambdas=min_lambdas,noise_mode=noise_mode)
        # else
            sol = simulate_steady_state(ode_func!, u0, tspan, params, saveat=20.0, n_tail = 100)
        # end
        push!(solutions, sol)
        
        # if i % 5 == 0
        #     print(".")
        # end
    end
    # println(" Done!")
    
    # Calculate thresholds
    if length(racipe_thresholds) == length(node_names)
        thresholds = copy(racipe_thresholds)
    else
        @warn "Improper threshodls. Cannot discretize. Exiting..."
        return
    end
    
    # Analyze each simulation
    # println("\nAnalyzing results...")
    for (i, (sol_matrix, t_vec)) in enumerate(solutions)
        # sol_matrix is nodes × tail_length, t_vec is the matching time vector
        # Discretize each timepoint
        discrete_states = [Tuple(sol_matrix[:, k] .> thresholds) for k in axes(sol_matrix, 2)]

        # Per-node stats over the tail
        stats = Dict(
            "Mean"     => [mean(row)     for row in eachrow(sol_matrix)],
            "StDev"    => [std(row)      for row in eachrow(sol_matrix)],
            "Skewness" => [skewness(row) for row in eachrow(sol_matrix)],
            "Kurtosis" => [kurtosis(row) for row in eachrow(sol_matrix)],
        )

        # Calculate metrics
        mrt = calculate_mrt(discrete_states)
        switches = count_switches(discrete_states)
        transitions = count_transitions(discrete_states)

        # Store result
        result = StochasticResult(
            discrete_states,
            thresholds,
            mrt,
            switches,
            transitions,
            sigma,
            stats
        )

        push!(results, result)
    end
    if !getSolutions
        return results
    else
        # Pre-allocate by building all rows at once
        all_times = Float64[]
        all_states = [Float64[] for _ in node_names]
        all_sim_ids = Int[]

        for (i, (states, t)) in enumerate(solutions)
            # states is nodes × tail_length, t is the matching time vector
            append!(all_times, t)
            append!(all_sim_ids, fill(i, length(t)))
            for (j, _) in enumerate(node_names)
                append!(all_states[j], states[j, :])
            end
        end

        sol_df = DataFrame(SimID = all_sim_ids, Time = all_times)
        for (j, node) in enumerate(node_names)
            sol_df[!, Symbol(node)] = all_states[j]
        end

        return results, sol_df
    end
end

"""
    aggregate_mrt(results)

Calculate mean MRT across multiple simulations

# Arguments
- `results`: Vector of StochasticResult objects

# Returns
- Dictionary of state => mean MRT
"""
function aggregate_mrt(results::Vector{StochasticResult})
    # Collect all unique states
    all_states = Set{Tuple}()
    for result in results
        union!(all_states, keys(result.mrt))
    end
    
    # Calculate mean for each state
    mean_mrt = Dict{Tuple, Float64}()
    
    for state in all_states
        mrts = [get(result.mrt, state, 0.0) for result in results]
        mean_mrt[state] = mean(mrts)
    end
    
    return mean_mrt
end

function aggregate_transitions(results::Vector{StochasticResult})
    all_states = Set{Tuple}()
    for result in results
        union!(all_states, keys(result.transition_counts))
    end
    
    # Calculate mean for each state
    mean_counts = Dict{Tuple, Float64}()
    
    for state in all_states
        counts = [get(result.transition_counts, state, 0.0) for result in results]
        mean_counts[state] = mean(counts)
    end
    return mean_counts
end

"""
    analyze_noise_effects(ode_func!, params, lambda_indices, node_names;
                         sigma_values=[0.0, 0.01, 0.05, 0.1, 0.2],
                         num_sims=25, tspan=(0.0, 1000.0))

Analyze effects of different noise levels

# Returns
- DataFrame with results for each noise level
"""
function analyze_noise_effects(ode_func!, params, param_names, lambda_indices, node_names;
                                u0List::Vector{Vector{Float64}} = Vector{Vector{Float64}}(),racipe_thresholds::Vector{Float64} = Vector{Float64}(),
                               sigma_values::Vector{Float64}=[0.0, 0.01, 0.05, 0.1, 0.2],
                               num_sims::Int=25,
                                            max_lambdas::Vector{Float64}=Float64[],
                                            min_lambdas::Vector{Float64}=Float64[],
                               tspan::Tuple{Float64,Float64}=(0.0, 1000.0), 
                               noise_mode::String = "Additive", saveat::Float64=1.0, 
                               getSolutions::Bool = false, 
                               dt_vals::Vector{Float64}=[0.01])
    
    all_results = Dict{Tuple{Float64, Float64}, Vector{StochasticResult}}()
    all_solutions = Dict{Tuple{Float64, Float64}, DataFrame}()
    # println("="^60)
    # println("ANALYZING NOISE EFFECTS")
    # println("="^60)
    # println("Noise levels: $sigma_values")
    # println("Simulations per level: $num_sims")
    # println("="^60)
    
    for sigma in sigma_values
        # println("\n--- sigma = $sigma ---")
        for dt in dt_vals
            if getSolutions
                results, solutions = run_multiple_stochastic_simulations(
                    ode_func!, params, param_names, lambda_indices, node_names,
                    num_sims=num_sims, tspan=tspan, sigma=sigma, dt=dt,
                    u0List = u0List, racipe_thresholds = racipe_thresholds,
                    max_lambdas=max_lambdas, min_lambdas=min_lambdas, noise_mode = noise_mode,
                    saveat=saveat, getSolutions=true
                )
                all_solutions[(sigma, dt)] = solutions
            else
                results = run_multiple_stochastic_simulations(
                    ode_func!, params, param_names, lambda_indices, node_names,
                    num_sims=num_sims, tspan=tspan, sigma=sigma, dt=dt,
                    u0List = u0List, racipe_thresholds = racipe_thresholds,
                    max_lambdas=max_lambdas, min_lambdas=min_lambdas, noise_mode = noise_mode,
                    saveat=saveat, getSolutions = false
                )
            end
            all_results[(sigma,dt)] = results
        end
        
    end
    # print(all_results)
    # Create summary DataFrame
    summary = create_summary_dataframe(all_results, node_names)
    GC.gc()
    if getSolutions
        return all_results, summary, all_solutions
    else
        return all_results, summary
    end
end

function analyze_noise_effects_det(ode_func!, params, param_names, lambda_indices, node_names;
                                u0List::Vector{Vector{Float64}} = Vector{Vector{Float64}}(),racipe_thresholds::Vector{Float64} = Vector{Float64}(),
                            #    sigma_values::Vector{Float64}=[0.0, 0.01, 0.05, 0.1, 0.2],                               
                               num_sims::Int=25,
                                            max_lambdas::Vector{Float64}=Float64[],
                                            min_lambdas::Vector{Float64}=Float64[],
                               tspan::Tuple{Float64,Float64}=(0.0, 1000.0), 
                               noise_mode::String = "Additive", saveat::Float64=1.0)
    
    # all_results = Dict{Float64, Vector{StochasticResult}}()
    
    # println("="^60)
    # println("ANALYZING NOISE EFFECTS")
    # println("="^60)
    # println("Noise levels: $sigma_values")
    # println("Simulations per level: $num_sims")
    # println("="^60)
    results = run_multiple_stochastic_simulations(
            ode_func!, params, param_names, lambda_indices, node_names,
            num_sims=num_sims, tspan=tspan, sigma=0.0,
            u0List = u0List, racipe_thresholds = racipe_thresholds,
            max_lambdas=max_lambdas, min_lambdas=min_lambdas, noise_mode = noise_mode,
            saveat=saveat
        )
    GC.gc()
    return results
    
    # for sigma in sigma_values
    #     # println("\n--- sigma = $sigma ---")
        
        
    #     all_results[sigma] = results
    # end
    # # print(all_results)
    # # Create summary DataFrame
    # summary = create_summary_dataframe(all_results, node_names)
    
    # return all_results, summary
end

"""
    analyze_noise_effects_parallel(ode_func!, params, param_names, lambda_indices, node_names;
                                   sigma_values=[0.0, 0.01, 0.05, 0.1, 0.2], ...)

Parallel version - distributes ALL simulations across threads
Parallelizes over (noise_level, simulation) combinations
"""
function analyze_noise_effects_parallel(
                                ode_func!, params, param_names, lambda_indices, node_names;
                                u0List::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
                                racipe_thresholds::Vector{Float64}=Vector{Float64}(),
                                sigma_values::Vector{Float64}=[0.0, 0.01, 0.05, 0.1, 0.2],
                                num_sims::Int=25,
                                max_lambdas::Vector{Float64}=Float64[],
                                min_lambdas::Vector{Float64}=Float64[],
                                tspan::Tuple{Float64,Float64}=(0.0, 1000.0),
                                noise_mode::String="Additive",
                                seed::Int=123, saveat::Float64=1.0, dt_vals::Vector{Float64}=[0.01])
    
    # Random.seed!(seed)
    n_nodes = length(node_names)
    n_noise_levels = length(sigma_values)
    total_sims = n_noise_levels * num_sims
    
    println("="^60)
    println("PARALLEL NOISE ANALYSIS")
    println("="^60)
    println("Noise levels: $sigma_values")
    println("Simulations per level: $num_sims")
    println("Total simulations: $total_sims")
    println("Threads available: $(Threads.nthreads())")
    println("="^60)
    
    # Prepare initial conditions (same for all noise levels)
    if isempty(u0List)
        randInit = true
    else
        compatibles = findall(x -> length(x) == length(node_names), u0List)
        if isempty(compatibles)
            randInit = true
        else
            u0List = u0List[compatibles]
            randInit = false
        end
    end
    
    if length(u0List) < num_sims
        u0List = u0List[mod1.(1:num_sims, length(u0List))]
    end
    
    # Generate all initial conditions
    all_u0 = Vector{Vector{Float64}}(undef, num_sims)
    for i in 1:num_sims
        if randInit
            u0 = zeros(n_nodes)
            for j in 1:n_nodes
                prod_idx = findfirst(x -> occursin("Prod_of_$(node_names[j])", x), param_names)
                deg_idx = findfirst(x -> occursin("Deg_of_$(node_names[j])", x), param_names)
                u0[j] = rand() * 1.5 * params[prod_idx] / params[deg_idx]
            end
            all_u0[i] = u0
        else
            all_u0[i] = u0List[i]
        end
    end
    
    # Create all (noise_level, sim_index, dt) combinations
    combinations = [(sigma, i, dt) for dt in dt_vals for sigma in sigma_values for i in 1:num_sims]
    n_combinations = length(combinations)
    
    println("\nRunning $n_combinations simulations in parallel...")
    
    # Preallocate storage for all results
    all_solutions = Vector{ODESolution}(undef, n_combinations)
    all_metadata = Vector{Tuple{Float64, Int}}(undef, n_combinations)  # (sigma, sim_index)
    
    # Parallel simulation over ALL combinations
    for idx in 1:n_combinations
        sigma, sim_idx, dt = combinations[idx]
        u0 = all_u0[sim_idx]
        
        all_solutions[idx] = simulate_with_noise(
            ode_func!, u0, tspan, params, lambda_indices;
            sigma=sigma, dt=dt, saveat=saveat,
            max_lambdas=max_lambdas, min_lambdas=min_lambdas, noise_mode=noise_mode
        )
        all_metadata[idx] = (sigma, sim_idx)
    end
    
    println("Done!")
    
    # Calculate thresholds
    if length(racipe_thresholds) == length(node_names)
        thresholds = copy(racipe_thresholds)
    else
        @warn "Improper thresholds. Cannot discretize."
        return Dict{Float64, Vector{StochasticResult}}(), DataFrame()
    end
    println("Converting solutions to arrays...")
    all_solution_arrays = Vector{Matrix{Float64}}(undef, n_combinations)

    for idx in 1:n_combinations
        sol = all_solutions[idx]
        # Convert once, serially
        all_solution_arrays[idx] = Array(sol)  # This is thread-safe
    end
    # Analyze results in parallel
    println("Analyzing results in parallel...")
    all_stochastic_results = Vector{StochasticResult}(undef, n_combinations)
    
    # Threads.@threads for idx in 1:n_combinations
    #     sol = all_solutions[idx]
    #     sigma, sim_idx = all_metadata[idx]
        
    #     discrete_states = discretize_trajectory(sol, thresholds)
    #     mrt = calculate_mrt(discrete_states)
    #     switches = count_switches(discrete_states)
        
    #     all_stochastic_results[idx] = StochasticResult(
    #         sol,
    #         collect(sol.t),
    #         states_matrix,
    #         discrete_states,
    #         thresholds,
    #         mrt,
    #         switches,
    #         sigma
    #     )
    # end

    Threads.@threads for idx in 1:n_combinations
        try
            sol = all_solutions[idx]
            sigma, sim_idx = all_metadata[idx]

            discrete_states = discretize_trajectory(sol, thresholds)
            n_cut = floor(Int, (1 - frac) * length(discrete_states))
            discrete_states = discrete_states[n_cut+1:end]
            mrt             = calculate_mrt(discrete_states)
            switches        = count_switches(discrete_states)
            transitions     = count_transitions(discrete_states)
            stats           = get_trajectory_stats(sol)

            all_stochastic_results[idx] = StochasticResult(
                discrete_states,
                thresholds,
                mrt,
                switches,
                transitions,
                sigma,
                stats
            )

        catch e
            @error "Analysis failed for simulation idx=$idx, param=$(all_metadata[idx])" exception=(e, catch_backtrace())
            rethrow(e)
        end
    end
    
    # Group results by noise level
    println("Grouping results by noise level...")
    results_by_noise = Dict{Float64, Vector{StochasticResult}}()
    
    for idx in 1:n_combinations
        sigma, sim_idx = all_metadata[idx]
        
        if !haskey(results_by_noise, sigma)
            results_by_noise[sigma] = StochasticResult[]
        end
        
        push!(results_by_noise[sigma], all_stochastic_results[idx])
    end
    
    # Create summary
    summary = create_summary_dataframe(results_by_noise, node_names)
    
    return results_by_noise, summary
end


using Distributed  # Add to top of module

"""
    analyze_noise_effects_distributed(...)

Uses Distributed.jl (pmap) instead of threading
Like R's mclapply - spawns separate processes for each simulation
"""
function analyze_noise_effects_distributed(
                                ode_func!, params, param_names, lambda_indices, node_names;
                                u0List::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
                                racipe_thresholds::Vector{Float64}=Vector{Float64}(),
                                sigma_values::Vector{Float64}=[0.0, 0.01, 0.05, 0.1, 0.2],
                                num_sims::Int=25,
                                min_lambdas::Vector{Float64}=Float64[],
                                max_lambdas::Vector{Float64}=Float64[],
                                tspan::Tuple{Float64,Float64}=(0.0, 1000.0),
                                dt::Float64=0.01,
                                saveat::Float64=1.0,
                                noise_mode::String="Additive",
                                seed::Int=123)
    
    # Random.seed!(seed)
    n_nodes = length(node_names)
    
    # ============================================================
    # PREPARE DATA (main process)
    # ============================================================
    
    println("Preparing initial conditions...")
    
    # Setup initial conditions
    use_racipe_ics = false
    if !isempty(u0List)
        compatibles = findall(x -> length(x) == length(node_names), u0List)
        if !isempty(compatibles)
            u0List = u0List[compatibles]
            use_racipe_ics = true
        end
    end
    
    if use_racipe_ics && length(u0List) < num_sims
        u0List = u0List[mod1.(1:num_sims, length(u0List))]
    end
    
    # Generate all initial conditions
    all_u0 = Vector{Vector{Float64}}(undef, num_sims)
    for i in 1:num_sims
        if use_racipe_ics
            all_u0[i] = u0List[i]
        else
            u0 = zeros(n_nodes)
            for j in 1:n_nodes
                prod_idx = findfirst(x -> occursin("Prod_of_$(node_names[j])", x), param_names)
                deg_idx = findfirst(x -> occursin("Deg_of_$(node_names[j])", x), param_names)
                u0[j] = rand() * 1.5 * params[prod_idx] / params[deg_idx]
            end
            all_u0[i] = u0
        end
    end
    
    # Validate thresholds
    if length(racipe_thresholds) != length(node_names)
        @warn "Improper thresholds. Cannot discretize."
        return Dict{Float64, Vector{StochasticResult}}(), DataFrame()
    end
    thresholds = copy(racipe_thresholds)
    
    # Create all (noise, sim) combinations
    combinations = [(sigma, i) for sigma in sigma_values for i in 1:num_sims]
    n_total = length(combinations)
    
    println("Distributing $n_total simulations across $(nprocs()-1) workers...")
    
    # ============================================================
    # DISTRIBUTED SIMULATION (pmap)
    # ============================================================
    
    # pmap automatically distributes across worker processes
    results_list = pmap(combinations, distributed=true, batch_size=1) do (sigma, sim_idx)
        u0 = all_u0[sim_idx]
        
        # Run simulation on worker
        sol = simulate_with_noise(
            ode_func!, u0, tspan, params, lambda_indices,
            sigma=sigma, dt=dt, saveat=saveat,
            max_lambdas=max_lambdas, min_lambdas=min_lambdas, noise_mode=noise_mode
        )
        
        # Discretize and analyze on worker
        discrete_states = discretize_trajectory(sol, thresholds)
        mrt         = calculate_mrt(discrete_states)
        switches    = count_switches(discrete_states)
        transitions = count_transitions(discrete_states)
        stats       = get_trajectory_stats(sol)

        # Return minimal data (avoid transferring huge ODESolution objects)
        return (
            sigma           = sigma,
            sim_idx         = sim_idx,
            mrt             = mrt,
            switches        = switches,
            transitions     = transitions,
            stats           = stats,
            discrete_states = discrete_states,
        )
    end

    println("All simulations complete. Assembling results...")

    # ============================================================
    # ASSEMBLE RESULTS (main process)
    # ============================================================

    # Group by noise level
    results_by_noise = Dict{Float64, Vector{StochasticResult}}()
    for sigma in sigma_values
        results_by_noise[sigma] = StochasticResult[]
    end

    for result in results_list
        stoch_result = StochasticResult(
            result.discrete_states,
            thresholds,
            result.mrt,
            result.switches,
            result.transitions,
            result.sigma,
            result.stats
        )

        push!(results_by_noise[result.sigma], stoch_result)
    end
    
    # Create summary
    summary = create_summary_dataframe(results_by_noise, node_names)
    
    return results_by_noise, summary
end

# Export it
export analyze_noise_effects_distributed

"""
    create_summary_dataframe(all_results, node_names)

Create summary DataFrame from analysis results
"""
function create_summary_dataframe(all_results::Dict{Tuple{Float64, Float64}, Vector{StochasticResult}},
                                 node_names::Vector{String})
    
    rows = []
    
    for (sigma, dt) in sort(collect(keys(all_results)))
        # Mean MRT
        results = all_results[(sigma, dt)]
        mean_mrt = aggregate_mrt(results)
        
        # Mean switches
        mean_switches = mean([r.switches for r in results])
        std_switches = std([r.switches for r in results])
        
        # Add rows for each state
        for (state, mrt) in sort(collect(mean_mrt), by=x->x[2], rev=true)
            push!(rows, (
                sigma = sigma,
                dt = dt,
                state = string(state),
                mrt = mrt,
                mean_switches = mean_switches,
                std_switches = std_switches
            ))
        end
    end
    
    return DataFrame(rows)
end
end # module