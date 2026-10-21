# Auto-generated ODE system from topology file
# Nodes: A, B, C, D

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
    D = u[4]

    # Extract parameters
    Prod_of_A = p[1]
    Prod_of_B = p[2]
    Prod_of_C = p[3]
    Prod_of_D = p[4]
    Deg_of_A = p[5]
    Deg_of_B = p[6]
    Deg_of_C = p[7]
    Deg_of_D = p[8]
    Trd_of_BToA = p[9]
    Num_of_BToA = p[10]
    Inh_of_BToA = p[11]
    Trd_of_CToA = p[12]
    Num_of_CToA = p[13]
    Inh_of_CToA = p[14]
    Trd_of_DToA = p[15]
    Num_of_DToA = p[16]
    Inh_of_DToA = p[17]
    Trd_of_AToA = p[18]
    Num_of_AToA = p[19]
    Act_of_AToA = p[20]
    Trd_of_AToB = p[21]
    Num_of_AToB = p[22]
    Inh_of_AToB = p[23]
    Trd_of_CToB = p[24]
    Num_of_CToB = p[25]
    Inh_of_CToB = p[26]
    Trd_of_DToB = p[27]
    Num_of_DToB = p[28]
    Inh_of_DToB = p[29]
    Trd_of_BToB = p[30]
    Num_of_BToB = p[31]
    Act_of_BToB = p[32]
    Trd_of_AToC = p[33]
    Num_of_AToC = p[34]
    Inh_of_AToC = p[35]
    Trd_of_BToC = p[36]
    Num_of_BToC = p[37]
    Inh_of_BToC = p[38]
    Trd_of_DToC = p[39]
    Num_of_DToC = p[40]
    Inh_of_DToC = p[41]
    Trd_of_CToC = p[42]
    Num_of_CToC = p[43]
    Act_of_CToC = p[44]
    Trd_of_AToD = p[45]
    Num_of_AToD = p[46]
    Inh_of_AToD = p[47]
    Trd_of_BToD = p[48]
    Num_of_BToD = p[49]
    Inh_of_BToD = p[50]
    Trd_of_CToD = p[51]
    Num_of_CToD = p[52]
    Inh_of_CToD = p[53]
    Trd_of_DToD = p[54]
    Num_of_DToD = p[55]
    Act_of_DToD = p[56]

    # ODEs for each node
    du[1] = Prod_of_A * H_shifted(B, Trd_of_BToA, Num_of_BToA, Inh_of_BToA) * H_shifted(C, Trd_of_CToA, Num_of_CToA, Inh_of_CToA) * H_shifted(D, Trd_of_DToA, Num_of_DToA, Inh_of_DToA) * H_shifted(A, Trd_of_AToA, Num_of_AToA, Act_of_AToA) - Deg_of_A * A
    du[2] = Prod_of_B * H_shifted(A, Trd_of_AToB, Num_of_AToB, Inh_of_AToB) * H_shifted(C, Trd_of_CToB, Num_of_CToB, Inh_of_CToB) * H_shifted(D, Trd_of_DToB, Num_of_DToB, Inh_of_DToB) * H_shifted(B, Trd_of_BToB, Num_of_BToB, Act_of_BToB) - Deg_of_B * B
    du[3] = Prod_of_C * H_shifted(A, Trd_of_AToC, Num_of_AToC, Inh_of_AToC) * H_shifted(B, Trd_of_BToC, Num_of_BToC, Inh_of_BToC) * H_shifted(D, Trd_of_DToC, Num_of_DToC, Inh_of_DToC) * H_shifted(C, Trd_of_CToC, Num_of_CToC, Act_of_CToC) - Deg_of_C * C
    du[4] = Prod_of_D * H_shifted(A, Trd_of_AToD, Num_of_AToD, Inh_of_AToD) * H_shifted(B, Trd_of_BToD, Num_of_BToD, Inh_of_BToD) * H_shifted(C, Trd_of_CToD, Num_of_CToD, Inh_of_CToD) * H_shifted(D, Trd_of_DToD, Num_of_DToD, Act_of_DToD) - Deg_of_D * D
end

function make_parameter_vector(param_dict::Dict)
    """
    Create parameter vector from dictionary
    param_dict should have keys matching parameter names
    """
    params = zeros(56)
    params[55] = get(param_dict, "Num_of_DToD", 0.0)
    params[21] = get(param_dict, "Trd_of_AToB", 0.0)
    params[10] = get(param_dict, "Num_of_BToA", 0.0)
    params[52] = get(param_dict, "Num_of_CToD", 0.0)
    params[20] = get(param_dict, "Act_of_AToA", 0.0)
    params[56] = get(param_dict, "Act_of_DToD", 0.0)
    params[11] = get(param_dict, "Inh_of_BToA", 0.0)
    params[53] = get(param_dict, "Inh_of_CToD", 0.0)
    params[44] = get(param_dict, "Act_of_CToC", 0.0)
    params[37] = get(param_dict, "Num_of_BToC", 0.0)
    params[14] = get(param_dict, "Inh_of_CToA", 0.0)
    params[29] = get(param_dict, "Inh_of_DToB", 0.0)
    params[34] = get(param_dict, "Num_of_AToC", 0.0)
    params[43] = get(param_dict, "Num_of_CToC", 0.0)
    params[42] = get(param_dict, "Trd_of_CToC", 0.0)
    params[49] = get(param_dict, "Num_of_BToD", 0.0)
    params[51] = get(param_dict, "Trd_of_CToD", 0.0)
    params[47] = get(param_dict, "Inh_of_AToD", 0.0)
    params[17] = get(param_dict, "Inh_of_DToA", 0.0)
    params[24] = get(param_dict, "Trd_of_CToB", 0.0)
    params[8] = get(param_dict, "Deg_of_D", 0.0)
    params[15] = get(param_dict, "Trd_of_DToA", 0.0)
    params[7] = get(param_dict, "Deg_of_C", 0.0)
    params[22] = get(param_dict, "Num_of_AToB", 0.0)
    params[30] = get(param_dict, "Trd_of_BToB", 0.0)
    params[13] = get(param_dict, "Num_of_CToA", 0.0)
    params[6] = get(param_dict, "Deg_of_B", 0.0)
    params[1] = get(param_dict, "Prod_of_A", 0.0)
    params[33] = get(param_dict, "Trd_of_AToC", 0.0)
    params[36] = get(param_dict, "Trd_of_BToC", 0.0)
    params[5] = get(param_dict, "Deg_of_A", 0.0)
    params[3] = get(param_dict, "Prod_of_C", 0.0)
    params[28] = get(param_dict, "Num_of_DToB", 0.0)
    params[48] = get(param_dict, "Trd_of_BToD", 0.0)
    params[16] = get(param_dict, "Num_of_DToA", 0.0)
    params[41] = get(param_dict, "Inh_of_DToC", 0.0)
    params[2] = get(param_dict, "Prod_of_B", 0.0)
    params[46] = get(param_dict, "Num_of_AToD", 0.0)
    params[26] = get(param_dict, "Inh_of_CToB", 0.0)
    params[27] = get(param_dict, "Trd_of_DToB", 0.0)
    params[9] = get(param_dict, "Trd_of_BToA", 0.0)
    params[25] = get(param_dict, "Num_of_CToB", 0.0)
    params[38] = get(param_dict, "Inh_of_BToC", 0.0)
    params[32] = get(param_dict, "Act_of_BToB", 0.0)
    params[4] = get(param_dict, "Prod_of_D", 0.0)
    params[35] = get(param_dict, "Inh_of_AToC", 0.0)
    params[50] = get(param_dict, "Inh_of_BToD", 0.0)
    params[54] = get(param_dict, "Trd_of_DToD", 0.0)
    params[45] = get(param_dict, "Trd_of_AToD", 0.0)
    params[12] = get(param_dict, "Trd_of_CToA", 0.0)
    params[19] = get(param_dict, "Num_of_AToA", 0.0)
    params[18] = get(param_dict, "Trd_of_AToA", 0.0)
    params[39] = get(param_dict, "Trd_of_DToC", 0.0)
    params[40] = get(param_dict, "Num_of_DToC", 0.0)
    params[23] = get(param_dict, "Inh_of_AToB", 0.0)
    params[31] = get(param_dict, "Num_of_BToB", 0.0)
    return params
end

function simulate_system(params::Vector{Float64}, u0::Vector{Float64}, tspan::Tuple)
    prob = ODEProblem(ode_system!, u0, tspan, params)
    sol = solve(prob, Tsit5(), saveat=0.1)
    return sol
end

# Node names for reference
NODE_NAMES = ["A", "B", "C", "D"]
N_NODES = 4
