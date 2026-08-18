"""
Generate ODE systems from RACIPE topology files
"""

module ODEgenerator

using DataFrames

export generate_ode_system, generate_ode_function

include("RACIPEdata.jl")
using .RACIPEdata

function get_all_nodes(topo::Topology)
    """Get unique list of all nodes in network"""
    nodes = unique(vcat(topo.sources, topo.targets))
    return sort(nodes)
end

function get_regulators(node::String, topo::Topology)
    """
    Get all regulators of a given node
    
    Returns: Vector of (source, type) tuples
    """
    regulators = []
    for i in 1:length(topo.targets)
        if topo.targets[i] == node
            push!(regulators, (topo.sources[i], topo.types[i]))
        end
    end
    return regulators
end

function generate_hill_functions()
    """
    Generate Hill function code
    Same as in models.jl
    """
    return """
# Hill function components
function H_minus(x, x0, n)
    return x0^n / (x0^n + x^n)
end

function H_plus(x, x0, n)
    return 1 - H_minus(x, x0, n)
end

function H_shifted(x, x0, n, lambda)
    if lambda>1 
        return H_minus(x, x0, n)/lambda + H_plus(x, x0, n)
    else
        return H_minus(x, x0, n) + lambda * H_plus(x, x0, n)
    end
end
"""
end

function generate_ode_system(topo::Topology, prs::ParameterRanges)
    """
    Generate ODE system code from topology
    
    Returns: String containing complete ODE function
    """
    nodes = get_all_nodes(topo)
    n_nodes = length(nodes)
    
    # Create node index mapping
    node_to_idx = Dict(node => i for (i, node) in enumerate(nodes))
    
    # Start building the function
    code = """
# Auto-generated ODE system from topology file
# Nodes: $(join(nodes, ", "))

using OrdinaryDiffEq
using DiffEqCallbacks

"""
    
    # Add Hill functions
    code *= generate_hill_functions()
    
    # Start ODE function
    code *= "\nfunction ode_system!(du, u, p, t)\n"
    code *= "    # Extract state variables\n"
    for (i, node) in enumerate(nodes)
        code *= "    $node = u[$i]\n"
    end
    code *= "\n"
    
    # Extract parameters
    code *= "    # Extract parameters\n"
    param_idx = 1
    param_mapping = Dict{String, Int}()
    
    for param_name in prs.names
        code *= "    $param_name = p[$param_idx]\n"
        param_mapping[param_name] = param_idx
        param_idx += 1
    end
    code *= "\n"
    
    # Generate ODE for each node
    code *= "    # ODEs for each node\n"
    for (idx, node) in enumerate(nodes)
        regulators = get_regulators(node, topo)
        
        prod_param = "Prod_of_$node"
        deg_param = "Deg_of_$node"
        
        # Start with production term
        if length(regulators) == 0
            # No regulators - constitutive production
            code *= "    du[$idx] = $prod_param - $deg_param * $node\n"
        else
            # Has regulators
            code *= "    du[$idx] = $prod_param"
            
            # Add regulatory terms
            for (regulator, reg_type) in regulators
                # Parameter names following RACIPE convention
                trd_param = "Trd_of_$(regulator)To$node"
                num_param = "Num_of_$(regulator)To$node"
                
                if reg_type == 2  # Inhibition
                    fold_param = "Inh_of_$(regulator)To$node"
                elseif reg_type == 1  # Activation
                    fold_param = "Act_of_$(regulator)To$node"
                else
                    error("Unknown regulation type: $reg_type")
                end
                
                code *= " * H_shifted($regulator, $trd_param, $num_param, $fold_param)"
            end
            
            # Add degradation
            code *= " - $deg_param * $node\n"
        end
    end
    
    code *= "end\n\n"
    
    # Add helper function to create parameter vector
    code *= "function make_parameter_vector(param_dict::Dict)\n"
    code *= "    \"\"\"\n"
    code *= "    Create parameter vector from dictionary\n"
    code *= "    param_dict should have keys matching parameter names\n"
    code *= "    \"\"\"\n"
    code *= "    params = zeros($(length(prs.names)))\n"
    for (param_name, idx) in param_mapping
        code *= "    params[$idx] = get(param_dict, \"$param_name\", 0.0)\n"
    end
    code *= "    return params\n"
    code *= "end\n\n"
    
    # Add helper to run simulation
    code *= """
function simulate_system(params::Vector{Float64}, u0::Vector{Float64}, tspan::Tuple)
    prob = ODEProblem(ode_system!, u0, tspan, params)
    sol = solve(prob, Tsit5(), saveat=0.1)
    return sol
end

# Node names for reference
NODE_NAMES = $(repr(nodes))
N_NODES = $n_nodes
"""
    
    return code
end

function generate_ode_function(topo_file::String, prs_file::String, 
                               output_file::String)
    """
    Generate ODE system and write to file
    
    Usage:
        generate_ode_function("TS.topo", "TS.prs", "TS_ode.jl")
    """
    topo = read_topo(topo_file)
    prs = read_prs(prs_file)
    
    code = generate_ode_system(topo, prs)
    
    open(output_file, "w") do f
        write(f, code)
    end
    
    println("ODE system written to $output_file")
    println("Nodes: $(get_all_nodes(topo))")
end

end # module
