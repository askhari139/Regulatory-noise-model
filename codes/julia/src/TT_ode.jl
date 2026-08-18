# Auto-generated ODE system from topology file
# Nodes: A, B, C

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
    C = u[3]

    # Extract parameters
    Prod_of_A = p[1]
    Prod_of_B = p[2]
    Prod_of_C = p[3]
    Deg_of_A = p[4]
    Deg_of_B = p[5]
    Deg_of_C = p[6]
    Trd_of_BToA = p[7]
    Num_of_BToA = p[8]
    Inh_of_BToA = p[9]
    Trd_of_CToA = p[10]
    Num_of_CToA = p[11]
    Inh_of_CToA = p[12]
    Trd_of_AToB = p[13]
    Num_of_AToB = p[14]
    Inh_of_AToB = p[15]
    Trd_of_CToB = p[16]
    Num_of_CToB = p[17]
    Inh_of_CToB = p[18]
    Trd_of_AToC = p[19]
    Num_of_AToC = p[20]
    Inh_of_AToC = p[21]
    Trd_of_BToC = p[22]
    Num_of_BToC = p[23]
    Inh_of_BToC = p[24]

    # ODEs for each node
    du[1] = Prod_of_A * H_shifted(B, Trd_of_BToA, Num_of_BToA, Inh_of_BToA) * H_shifted(C, Trd_of_CToA, Num_of_CToA, Inh_of_CToA) - Deg_of_A * A
    du[2] = Prod_of_B * H_shifted(A, Trd_of_AToB, Num_of_AToB, Inh_of_AToB) * H_shifted(C, Trd_of_CToB, Num_of_CToB, Inh_of_CToB) - Deg_of_B * B
    du[3] = Prod_of_C * H_shifted(A, Trd_of_AToC, Num_of_AToC, Inh_of_AToC) * H_shifted(B, Trd_of_BToC, Num_of_BToC, Inh_of_BToC) - Deg_of_C * C
end

function make_parameter_vector(param_dict::Dict)
    """
    Create parameter vector from dictionary
    param_dict should have keys matching parameter names
    """
    params = zeros(24)
    params[5] = get(param_dict, "Deg_of_B", 0.0)
    params[24] = get(param_dict, "Inh_of_BToC", 0.0)
    params[1] = get(param_dict, "Prod_of_A", 0.0)
    params[13] = get(param_dict, "Trd_of_AToB", 0.0)
    params[21] = get(param_dict, "Inh_of_AToC", 0.0)
    params[8] = get(param_dict, "Num_of_BToA", 0.0)
    params[16] = get(param_dict, "Trd_of_CToB", 0.0)
    params[19] = get(param_dict, "Trd_of_AToC", 0.0)
    params[22] = get(param_dict, "Trd_of_BToC", 0.0)
    params[4] = get(param_dict, "Deg_of_A", 0.0)
    params[9] = get(param_dict, "Inh_of_BToA", 0.0)
    params[3] = get(param_dict, "Prod_of_C", 0.0)
    params[10] = get(param_dict, "Trd_of_CToA", 0.0)
    params[2] = get(param_dict, "Prod_of_B", 0.0)
    params[23] = get(param_dict, "Num_of_BToC", 0.0)
    params[6] = get(param_dict, "Deg_of_C", 0.0)
    params[12] = get(param_dict, "Inh_of_CToA", 0.0)
    params[15] = get(param_dict, "Inh_of_AToB", 0.0)
    params[14] = get(param_dict, "Num_of_AToB", 0.0)
    params[18] = get(param_dict, "Inh_of_CToB", 0.0)
    params[20] = get(param_dict, "Num_of_AToC", 0.0)
    params[7] = get(param_dict, "Trd_of_BToA", 0.0)
    params[17] = get(param_dict, "Num_of_CToB", 0.0)
    params[11] = get(param_dict, "Num_of_CToA", 0.0)
    return params
end

function simulate_system(params::Vector{Float64}, u0::Vector{Float64}, tspan::Tuple)
    prob = ODEProblem(ode_system!, u0, tspan, params)
    sol = solve(prob, Tsit5(), saveat=0.1)
    return sol
end

# Node names for reference
NODE_NAMES = ["A", "B", "C"]
N_NODES = 3
