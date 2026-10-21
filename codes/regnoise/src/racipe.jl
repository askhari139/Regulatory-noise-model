# RACIPE inputs and the ODE model (SPEC §2, §5).

struct Edge
    src::Int          # source node index
    tgt::Int          # target node index
    act::Bool         # true = activation, false = inhibition
    name::String      # e.g. "BToA"
end

"""
Network model. `param_names` follows the .prs order, which is also the column order of
`<NET>_parameters.dat`, so a parameter vector is a row of that file.
"""
struct Network
    name::String
    nodes::Vector{String}
    edges::Vector{Edge}
    param_names::Vector{String}
    g_idx::Vector{Int}        # Prod_of_<node>
    k_idx::Vector{Int}        # Deg_of_<node>
    trd_idx::Vector{Int}      # per edge
    num_idx::Vector{Int}
    lam_idx::Vector{Int}
    in_ptr::Vector{Int}       # CSR: incoming edges of node i are in_edges[in_ptr[i]:in_ptr[i+1]-1]
    in_edges::Vector{Int}
end

nnodes(net::Network) = length(net.nodes)
nedges(net::Network) = length(net.edges)

function read_prs(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = split(lines[1], '\t')
    header[1] == "Parameter" || error("$path: unexpected header $(lines[1])")
    names = String[]; lo = Float64[]; hi = Float64[]
    for l in lines[2:end]
        f = split(l, '\t')
        length(f) == 4 || error("$path: malformed line '$l'")
        push!(names, f[1]); push!(lo, parse(Float64, f[2])); push!(hi, parse(Float64, f[3]))
    end
    return names, lo, hi
end

function read_topo(path)
    lines = filter(!isempty, strip.(readlines(path)))
    split(lines[1]) == ["Source", "Target", "Type"] || error("$path: unexpected header $(lines[1])")
    edges = Tuple{String,String,Int}[]
    for l in lines[2:end]
        f = split(l)
        length(f) == 3 || error("$path: malformed line '$l'")
        t = parse(Int, f[3])
        t in (1, 2) || error("$path: edge type must be 1 or 2, got $t")
        push!(edges, (f[1], f[2], t))
    end
    return edges
end

function load_network(cfg::Config, netname::AbstractString)
    names, _, _ = read_prs(joinpath(cfg.data_dir, netname * ".prs"))
    topo = read_topo(joinpath(cfg.data_dir, netname * ".topo"))
    pidx(n) = (i = findfirst(==(n), names); i === nothing ? error("$netname.prs lacks parameter $n") : i)

    nodes = [replace(n, "Prod_of_" => "") for n in names if startswith(n, "Prod_of_")]
    nid(n) = (i = findfirst(==(n), nodes); i === nothing ? error("$netname.topo: unknown node $n") : i)

    edges = Edge[]; trd = Int[]; num = Int[]; lam = Int[]
    for (s, t, ty) in topo
        en = "$(s)To$(t)"
        act = ty == 1
        push!(edges, Edge(nid(s), nid(t), act, en))
        push!(trd, pidx("Trd_of_$en")); push!(num, pidx("Num_of_$en"))
        push!(lam, pidx((act ? "Act_of_" : "Inh_of_") * en))
    end
    # every fold-change parameter in the .prs must belong to a topology edge
    nfold = count(n -> startswith(n, "Act_of_") || startswith(n, "Inh_of_"), names)
    nfold == length(edges) || error("$netname: $nfold fold-change parameters but $(length(edges)) edges")

    order = sortperm([e.tgt for e in edges])
    in_ptr = zeros(Int, length(nodes) + 1)
    in_ptr[1] = 1
    for e in edges
        in_ptr[e.tgt+1] += 1
    end
    cumsum!(in_ptr, in_ptr)
    return Network(String(netname), nodes, edges, names,
        [pidx("Prod_of_$n") for n in nodes], [pidx("Deg_of_$n") for n in nodes],
        trd, num, lam, in_ptr, order)
end

"Callable ODE right-hand side; λ values are read from `p`, so callbacks can change them."
struct RHS
    net::Network
end

function (f::RHS)(du, u, p, t)
    net = f.net
    @inbounds for i in 1:length(net.nodes)
        prod = 1.0
        for j in net.in_ptr[i]:(net.in_ptr[i+1]-1)
            e = net.in_edges[j]
            x = u[net.edges[e].src]
            x0n = p[net.trd_idx[e]]^p[net.num_idx[e]]
            hm = x0n / (x0n + x^p[net.num_idx[e]])
            λ = p[net.lam_idx[e]]
            prod *= net.edges[e].act ? (hm + λ * (1 - hm)) / λ : hm + λ * (1 - hm)
        end
        du[i] = p[net.g_idx[i]] * prod - p[net.k_idx[i]] * u[i]
    end
    return nothing
end

"Per-edge λ bounds from the config."
lambda_bounds(cfg::Config, net::Network) =
    ([e.act ? cfg.act_bounds[1] : cfg.inh_bounds[1] for e in net.edges],
     [e.act ? cfg.act_bounds[2] : cfg.inh_bounds[2] for e in net.edges])

"Read `<NET>_parameters.dat` → Dict ParamID => parameter vector (.prs order)."
function read_parameters(cfg::Config, net::Network)
    path = joinpath(cfg.data_dir, net.name * "_parameters.dat")
    np = length(net.param_names)
    lo, hi = lambda_bounds(cfg, net)
    out = Dict{Int,Vector{Float64}}()
    for l in eachline(path)
        isempty(strip(l)) && continue
        f = split(l, '\t')
        length(f) == np + 2 || error("$path: expected $(np + 2) columns, got $(length(f))")
        id = parse(Int, f[1])
        p = parse.(Float64, f[3:end])
        for (k, li) in enumerate(net.lam_idx)
            lo[k] - 1e-9 <= p[li] <= hi[k] + 1e-9 ||
                error("$path ParamID $id: λ $(net.param_names[li]) = $(p[li]) outside configured bounds [$(lo[k]), $(hi[k])]")
        end
        out[id] = p
    end
    return out
end

"Read `<NET>_solution.dat` → rows (ParamID, frequency, log2 expression vector)."
function read_solutions(cfg::Config, net::Network)
    path = joinpath(cfg.data_dir, net.name * "_solution.dat")
    n = nnodes(net)
    rows = Tuple{Int,Float64,Vector{Float64}}[]
    for l in eachline(path)
        isempty(strip(l)) && continue
        f = split(l, '\t')
        length(f) == n + 3 || error("$path: expected $(n + 3) columns, got $(length(f))")
        push!(rows, (parse(Int, f[1]), parse(Float64, f[3]), parse.(Float64, f[4:end])))
    end
    return rows
end

"θ_i = 2^(frequency-weighted mean of log2 expression over all steady states) (SPEC §5)."
function thresholds(rows, n::Int)
    w = sum(r[2] for r in rows)
    return [2.0^(sum(r[2] * r[3][i] for r in rows) / w) for i in 1:n]
end

"State code: bit i-1 set when node i is above threshold."
@inline function state_code(x, θ, m=1.0)
    c = 0
    @inbounds for i in eachindex(θ)
        x[i] > m * θ[i] && (c |= 1 << (i - 1))
    end
    return c
end

"Human-readable state label in node order, e.g. A high, B low → \"10\"."
state_label(code::Int, n::Int) = join(((code >> (i - 1)) & 1 for i in 1:n))
