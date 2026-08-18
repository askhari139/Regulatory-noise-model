# Auto-generated ODE system from topology file
# Nodes: A, B

using OrdinaryDiffEq
using DiffEqCallbacks

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

function ode_system!(du, u, p, t)
    # Extract state variables
    A = u[1]
    B = u[2]

    # Extract parameters
    Prod_of_A = p[1]
    Prod_of_B = p[2]
    Deg_of_A = p[3]
    Deg_of_B = p[4]
    Trd_of_AToB = p[5]
    Num_of_AToB = p[6]
    Inh_of_AToB = p[7]

    # ODEs for each node
    du[1] = Prod_of_A - Deg_of_A * A
    du[2] = Prod_of_B * H_shifted(A, Trd_of_AToB, Num_of_AToB, Inh_of_AToB) - Deg_of_B * B
end

function make_parameter_vector(param_dict::Dict)
    """
    Create parameter vector from dictionary
    param_dict should have keys matching parameter names
    """
    params = zeros(7)
    params[2] = get(param_dict, "Prod_of_B", 0.0)
    params[1] = get(param_dict, "Prod_of_A", 0.0)
    params[3] = get(param_dict, "Deg_of_A", 0.0)
    params[5] = get(param_dict, "Trd_of_AToB", 0.0)
    params[7] = get(param_dict, "Inh_of_AToB", 0.0)
    params[6] = get(param_dict, "Num_of_AToB", 0.0)
    params[4] = get(param_dict, "Deg_of_B", 0.0)
    return params
end

function simulate_system(params::Vector{Float64}, u0::Vector{Float64}, tspan::Tuple)
    prob = ODEProblem(ode_system!, u0, tspan, params)
    sol = solve(prob, Tsit5(), saveat=0.1)
    return sol
end

# Node names for reference
NODE_NAMES = ["A", "B"]
N_NODES = 2
