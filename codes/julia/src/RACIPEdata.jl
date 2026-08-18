"""
Core data structures and parsers for RACIPE files
"""

module RACIPEdata

using DataFrames
using CSV
using Statistics
using StatsBase

export ParameterRanges, Topology, Solutions, Parameters
export read_prs, read_topo, read_solutions, read_parameters
export get_mean_expression, filter_parameters
export discretize_racipe_states, identify_attractors, classify_attractors_by_stability
export sample_balanced_parameters, get_parameter_attractor_info
# Data structures
struct ParameterRanges
    """Structure to hold parameter ranges from .prs file"""
    names::Vector{String}
    min_values::Vector{Float64}
    max_values::Vector{Float64}
    reg_types::Vector{Int}
    
    function ParameterRanges(names, mins, maxs, regs)
        @assert length(names) == length(mins) == length(maxs) == length(regs)
        new(names, mins, maxs, regs)
    end
end

struct Topology
    """Structure to hold network topology from .topo file"""
    sources::Vector{String}
    targets::Vector{String}
    types::Vector{Int}  # 1 = activation, 2 = inhibition
    
    function Topology(sources, targets, types)
        @assert length(sources) == length(targets) == length(types)
        new(sources, targets, types)
    end
end

struct Solutions
    """Structure to hold RACIPE solutions from _solution.dat file"""
    data::DataFrame
    node_names::Vector{String}
    
    function Solutions(data, node_names)
        new(data, node_names)
    end
end

struct Parameters
    """Structure to hold parameter values from _parameters.dat file"""
    data::DataFrame
    param_names::Vector{String}
    
    function Parameters(data, param_names)
        new(data, param_names)
    end
end

# Parsing functions
function read_prs(filename::String)
    """
    Read parameter ranges file (.prs)
    
    Returns: ParameterRanges object
    """
    df = CSV.read(filename, DataFrame, delim='\t', comment="#")
    
    names = df.Parameter
    mins = df.Minimum_value
    maxs = df.Maximum_Value
    regs = df.Regulation_type
    
    return ParameterRanges(names, mins, maxs, regs)
end

function read_topo(filename::String)
    """
    Read topology file (.topo)
    
    Returns: Topology object
    """
    df = CSV.read(filename, DataFrame, delim=' ', comment="#")
    
    sources = string.(df.Source)
    targets = string.(df.Target)
    types = df.Type
    
    return Topology(sources, targets, types)
end

function extract_node_names(prs::ParameterRanges)
    """
    Extract node names from parameter names
    E.g., "Prod_of_A" -> "A"
    """
    nodes = Vector{String}()
    
    for param in prs.names
        if startswith(param, "Prod_of_")
            node = split(param, "_of_")[2]
            push!(nodes, node)
        end
    end
    
    return nodes
end

function read_solutions(filename::String, prs::ParameterRanges)
    """
    Read solutions file (_solution.dat)
    
    Columns:
    1: Parameter ID
    2: State number
    3: Frequency (% of initial conditions)
    4+: Log2 expression values for each node
    
    Returns: Solutions object
    """
    # Read without header
    df = CSV.read(filename, DataFrame, delim='\t', header=false)
    
    # Get node names
    node_names = extract_node_names(prs)
    n_nodes = length(node_names)
    
    # Rename columns
    col_names = ["ParamID", "NumStates", "Frequency"]
    append!(col_names, node_names)
    
    # Check if we have the right number of columns
    expected_cols = 3 + n_nodes
    if ncol(df) != expected_cols
        error("Expected $expected_cols columns but got $(ncol(df))")
    end
    
    rename!(df, Symbol.(col_names))
    
    return Solutions(df, node_names)
end

function read_parameters(filename::String, prs::ParameterRanges)
    """
    Read parameters file (_parameters.dat)
    
    Columns:
    1: Parameter ID
    2: Number of steady states
    3+: Parameter values in order from .prs file
    
    Returns: Parameters object
    """
    # Read without header
    df = CSV.read(filename, DataFrame, delim='\t', header=false)
    
    # Create column names
    col_names = ["ParamID", "NumStates"]
    append!(col_names, prs.names)
    
    # Check column count
    expected_cols = 2 + length(prs.names)
    if ncol(df) != expected_cols
        error("Expected $expected_cols columns but got $(ncol(df))")
    end
    
    rename!(df, Symbol.(col_names))
    
    return Parameters(df, prs.names)
end

# Analysis functions
function get_mean_expression(solutions::Solutions; 
                            weight_by_frequency::Bool=true, log_mean::Bool = true)
    """
    Calculate mean expression of all nodes across all parameters
    
    Parameters:
    - solutions: Solutions object
    - weight_by_frequency: If true, weight by convergence frequency
    
    Returns: DataFrame with mean expression for each node
    """
    node_names = solutions.node_names
    
    if weight_by_frequency
        # Weight each state by its frequency
        means = Dict{String, Float64}()
        
        for node in node_names
            # Calculate weighted mean
            # Group by parameter ID and calculate weighted average
            weighted_sum = 0.0
            total_weight = 0.0
            
            for row in eachrow(solutions.data)
                weight = row.Frequency / 100.0  # Convert percentage to fraction
                if !log_mean
                    value = (2).^row[Symbol(node)]
                else
                    value = row[Symbol(node)]
                end
                weighted_sum += weight * value
                total_weight += weight
            end
            
            means[node] = weighted_sum / total_weight
        end
        
    else
        # Simple average across all states
        means = Dict{String, Float64}()
        
        for node in node_names
            means[node] = mean(solutions.data[!, Symbol(node)])
        end
    end
    
    # Create DataFrame
    result = DataFrame(
        Node = collect(keys(means)),
        MeanExpression = [2^x for x in collect(values(means))]
    )
    # par_ids = sort(unique(solutions.data.ParamID))
    # df = DataFrame(ParamID = par_ids)
    # for nd in solutions.node_names
    #     df[!, Symbol(nd)] .= first(result[result.Node == nd, :MeanExpression])
    # end
    
    # return df
    return result
end

function get_mean_expression(parameters::Parameters)
    par_names = parameters.param_names
    par_ids = parameters.data.ParamID
    nodes = filter(n -> startswith(n, "Prod_of_"), par_names)
    nodes = replace.(nodes, "Prod_of_" => "")
    df = DataFrame(ParamID = par_ids)
    for node in nodes
        thresholds = filter(n -> startswith(n, "Trd_of_$(node)"), par_names)
        thresholds = parameters.data[!, Symbol.(thresholds)]
        mean_expressions = [sum(x)/length(x) for x in eachrow(thresholds)]
        df[!, Symbol(node)] = mean_expressions
    end
    return df
end
function get_mean_expression_per_parameter(solutions::Solutions; log_mean::Bool = true)
    """
    Calculate mean expression for each parameter set separately
    
    Returns: DataFrame with ParamID and mean expression for each node
    """
    node_names = solutions.node_names
    
    # Group by parameter ID
    grouped = groupby(solutions.data, :ParamID)
    
    results = DataFrame()
    
    for group in grouped
        param_id = first(group.ParamID)
        
        row_data = Dict("ParamID" => param_id)
        
        for node in node_names
            # Weighted mean within this parameter set
            if log_mean
                weighted_sum = sum(group.Frequency .* group[!, Symbol(node)])
            else
                weighted_sum = sum(group.Frequency .* (2).^group[!, Symbol(node)])
            end
            total_freq = sum(group.Frequency)
            row_data[node] = weighted_sum / total_freq
        end
        
        push!(results, row_data, cols=:union)
    end
    
    return results
end

function filter_parameters(params::Parameters, solutions::Solutions;
                          num_states::Union{Int,Nothing}=nothing,
                          num_states_range::Union{Tuple{Int,Int},Nothing}=nothing,
                          prod_deg_ratio::Union{Dict{String,Float64},Nothing}=nothing,
                          custom_filter::Union{Function,Nothing}=nothing)
    """
    Filter parameter sets based on various criteria
    
    Parameters:
    - params: Parameters object
    - solutions: Solutions object
    - num_states: Exact number of steady states (e.g., 2 for bistable)
    - num_states_range: Range of steady states (min, max)
    - prod_deg_ratio: Dict of node => minimum ratio of production/degradation
    - custom_filter: Custom function that takes a row and returns Bool
    
    Returns: Filtered Parameters and Solutions objects
    """
    # Start with all parameter IDs
    mask = trues(nrow(params.data))
    
    # Filter by number of states
    if num_states !== nothing
        mask .&= params.data.NumStates .== num_states
    end
    
    if num_states_range !== nothing
        min_states, max_states = num_states_range
        mask .&= (params.data.NumStates .>= min_states) .& 
                 (params.data.NumStates .<= max_states)
    end
    
    # Filter by production/degradation ratio
    if prod_deg_ratio !== nothing
        for (node, min_ratio) in prod_deg_ratio
            prod_col = Symbol("Prod_of_$node")
            deg_col = Symbol("Deg_of_$node")
            
            if hasproperty(params.data, prod_col) && hasproperty(params.data, deg_col)
                ratio = params.data[!, prod_col] ./ params.data[!, deg_col]
                mask .&= ratio .>= min_ratio
            else
                @warn "Node $node not found in parameters"
            end
        end
    end
    
    # Apply custom filter
    if custom_filter !== nothing
        for (i, row) in enumerate(eachrow(params.data))
            if mask[i]  # Only check if not already filtered out
                mask[i] = custom_filter(row)
            end
        end
    end
    
    # Filter parameter data
    filtered_params_data = params.data[mask, :]
    filtered_params = Parameters(filtered_params_data, params.param_names)
    
    # Filter solutions data (keep only solutions from filtered parameters)
    filtered_param_ids = Set(filtered_params_data.ParamID)
    sol_mask = [id in filtered_param_ids for id in solutions.data.ParamID]
    filtered_sol_data = solutions.data[sol_mask, :]
    filtered_solutions = Solutions(filtered_sol_data, solutions.node_names)
    
    println("Filtered: $(sum(mask)) / $(length(mask)) parameter sets remain")
    
    return filtered_params, filtered_solutions
end

"""
    discretize_racipe_states(solutions::Solutions, thresholds::Vector{Float64})

Discretize RACIPE log2 expression values to binary states

# Arguments
- solutions: Solutions object
- thresholds: Vector of thresholds (one per node), typically from mean expression

# Returns
- DataFrame with ParamID, StateNum, Frequency, and DiscreteState columns
"""
function discretize_racipe_states(solutions::Solutions; thresholds::Vector{Float64} = nothing, 
    log_mean::Bool = true)
    
    @assert length(thresholds) == length(solutions.node_names)
    
    if isnothing(thresholds)
        thresholds = get_mean_expression(solutions; log_mean = log_mean)
    end
    # Add discrete state column
    discrete_df = copy(solutions.data)
    
    # Convert log2 expression to binary
    discrete_states = []
    for row in eachrow(discrete_df)
        state = Tuple(2^row[Symbol(node)] > thresholds[i] ? 1 : 0 
                     for (i, node) in enumerate(solutions.node_names))
        push!(discrete_states, string(state))
    end
    
    discrete_df.DiscreteState = discrete_states
    
    return discrete_df
end

"""
    identify_attractors(solutions::Solutions, thresholds::Vector{Float64})

Identify attractor structures and group parameters

# Returns
- Dictionary mapping attractor_signature => attractor_info
  where attractor_signature is Set of discrete states (e.g., Set(["(0,1)", "(1,0)"]))
  and attractor_info contains:
    - param_ids: Vector of parameter IDs with this attractor
    - states: Vector of states in this attractor
    - frequency: Total frequency across all parameters
    - balance_scores: Dict of param_id => balance score
"""
function identify_attractors(solutions::Solutions; thresholds::Vector{Float64}=nothing)
    
    # Discretize first
    discrete_df = discretize_racipe_states(solutions; thresholds = thresholds)
    
    # Group by parameter
    attractors = Dict{Set{String}, Dict}()
    
    for param_id in unique(discrete_df.ParamID)
        param_data = discrete_df[discrete_df.ParamID .== param_id, :]
        
        # Get all states for this parameter
        states = Set(param_data.DiscreteState)
        
        # Calculate balance score for this parameter
        # For bistable: 1.0 = perfectly balanced, 0.0 = one state dominates
        # For monostable: always 1.0
        # For tristable: balance of least vs most common
        freqs = param_data.Frequency
        balance = length(states) == 1 ? 1.0 : minimum(freqs) / maximum(freqs)
        
        # Total frequency
        total_freq = sum(freqs)
        
        # Create or update attractor entry
        if !haskey(attractors, states)
            attractors[states] = Dict(
                "param_ids" => Int[],
                "states" => collect(states),
                "frequencies" => Float64[],
                "balance_scores" => Dict{Int, Float64}(),
                "state_distributions" => Dict{Int, Dict{String, Float64}}()
            )
        end
        
        push!(attractors[states]["param_ids"], param_id)
        push!(attractors[states]["frequencies"], total_freq)
        attractors[states]["balance_scores"][param_id] = balance
        
        # Store full state distribution for this parameter
        state_dist = Dict{String, Float64}()
        for row in eachrow(param_data)
            state_dist[row.DiscreteState] = row.Frequency
        end
        attractors[states]["state_distributions"][param_id] = state_dist
    end
    
    return attractors
end

"""
    classify_attractors_by_stability(attractors::Dict)

Classify attractors by stability type (mono, bi, tri, tetra, etc.)

# Returns
- Dictionary mapping stability_class => attractor_info
"""
function classify_attractors_by_stability(attractors::Dict)
    
    classified = Dict{String, Vector}()
    
    for (states, info) in attractors
        n_states = length(states)
        
        class = if n_states == 1
            "monostable"
        elseif n_states == 2
            "bistable"
        elseif n_states == 3
            "tristable"
        elseif n_states == 4
            "tetrastable"
        else
            "multistable_$(n_states)"
        end
        
        if !haskey(classified, class)
            classified[class] = []
        end
        
        push!(classified[class], (states=states, info=info))
    end
    
    return classified
end

"""
    sample_balanced_parameters(attractors::Dict, stability_class::String, 
                               n_samples::Int; 
                               weight_by_frequency::Bool=true,
                               prefer_balanced::Bool=true)

Sample parameters from an attractor class, preferring balanced ones

# Arguments
- attractors: Output from identify_attractors
- stability_class: "monostable", "bistable", etc.
- n_samples: Number of parameters to sample
- weight_by_frequency: Weight by attractor frequency
- prefer_balanced: Prioritize balanced parameters (equal state frequencies)

# Returns
- Vector of parameter IDs
"""
function sample_balanced_parameters(attractors::Dict, 
                                   stability_class::String,
                                   n_samples::Int;
                                   weight_by_frequency::Bool=true,
                                   prefer_balanced::Bool=true)
    
    # Classify attractors
    classified = classify_attractors_by_stability(attractors)
    
    if !haskey(classified, stability_class)
        @warn "No $(stability_class) attractors found"
        return Int[]
    end
    
    # Collect all parameters from this class
    all_param_ids = Int[]
    all_frequencies = Float64[]
    all_balance_scores = Float64[]
    
    for (states, info) in classified[stability_class]
        append!(all_param_ids, info["param_ids"])
        append!(all_frequencies, info["frequencies"])
        
        for param_id in info["param_ids"]
            push!(all_balance_scores, info["balance_scores"][param_id])
        end
    end
    
    if isempty(all_param_ids)
        return Int[]
    end
    
    # Calculate sampling weights
    if prefer_balanced && weight_by_frequency
        # Combine balance score and frequency
        weights = all_balance_scores .* (all_frequencies ./ sum(all_frequencies))
    elseif prefer_balanced
        weights = all_balance_scores
    elseif weight_by_frequency
        weights = all_frequencies
    else
        weights = ones(length(all_param_ids))
    end
    
    # Normalize weights
    weights = weights ./ sum(weights)
    
    # Sample with replacement based on weights
    n_available = length(all_param_ids)
    n_to_sample = min(n_samples, n_available)
    
    # Sample without replacement, weighted
    sampled_indices = sample(1:n_available, Weights(weights), 
                            n_to_sample, replace=false)
    
    return all_param_ids[sampled_indices]
end

"""
    get_parameter_attractor_info(param_id::Int, attractors::Dict)

Get attractor information for a specific parameter

# Returns
- Dictionary with states, balance_score, frequency
"""
function get_parameter_attractor_info(param_id::Int, attractors::Dict)
    
    for (states, info) in attractors
        if param_id in info["param_ids"]
            idx = findfirst(==(param_id), info["param_ids"])
            
            return Dict(
                "states" => collect(states),
                "balance_score" => info["balance_scores"][param_id],
                "frequency" => info["frequencies"][idx],
                "state_distribution" => info["state_distributions"][param_id]
            )
        end
    end
    
    return nothing
end

# Export new functions


end # module
