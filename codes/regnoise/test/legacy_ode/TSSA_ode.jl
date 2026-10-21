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
    Trd_of_BToA = p[5]
    Num_of_BToA = p[6]
    Inh_of_BToA = p[7]
    Trd_of_AToA = p[8]
    Num_of_AToA = p[9]
    Act_of_AToA = p[10]
    Trd_of_AToB = p[11]
    Num_of_AToB = p[12]
    Inh_of_AToB = p[13]
    Trd_of_BToB = p[14]
    Num_of_BToB = p[15]
    Act_of_BToB = p[16]

    # ODEs for each node
    du[1] = Prod_of_A * H_shifted(B, Trd_of_BToA, Num_of_BToA, Inh_of_BToA) * H_shifted(A, Trd_of_AToA, Num_of_AToA, Act_of_AToA) - Deg_of_A * A
    du[2] = Prod_of_B * H_shifted(A, Trd_of_AToB, Num_of_AToB, Inh_of_AToB) * H_shifted(B, Trd_of_BToB, Num_of_BToB, Act_of_BToB) - Deg_of_B * B
end

function make_parameter_vector(param_dict::Dict)
    """
    Create parameter vector from dictionary
    param_dict should have keys matching parameter names
    """
    params = zeros(16)
    params[1] = get(param_dict, "Prod_of_A", 0.0)
    params[11] = get(param_dict, "Trd_of_AToB", 0.0)
    params[16] = get(param_dict, "Act_of_BToB", 0.0)
    params[6] = get(param_dict, "Num_of_BToA", 0.0)
    params[3] = get(param_dict, "Deg_of_A", 0.0)
    params[10] = get(param_dict, "Act_of_AToA", 0.0)
    params[7] = get(param_dict, "Inh_of_BToA", 0.0)
    params[9] = get(param_dict, "Num_of_AToA", 0.0)
    params[2] = get(param_dict, "Prod_of_B", 0.0)
    params[8] = get(param_dict, "Trd_of_AToA", 0.0)
    params[13] = get(param_dict, "Inh_of_AToB", 0.0)
    params[12] = get(param_dict, "Num_of_AToB", 0.0)
    params[14] = get(param_dict, "Trd_of_BToB", 0.0)
    params[15] = get(param_dict, "Num_of_BToB", 0.0)
    params[5] = get(param_dict, "Trd_of_BToA", 0.0)
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
