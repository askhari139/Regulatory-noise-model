# Parameter-set selection (SPEC §4). Run once per network; every mode reuses it.

const CLASS_NAMES = Dict(1 => "monostable", 2 => "bistable", 3 => "tristable", 4 => "tetrastable")

selection_path(cfg::Config, net::AbstractString) = joinpath(run_dir(cfg), "selection", net * ".csv")
thresholds_path(cfg::Config, net::AbstractString) = joinpath(run_dir(cfg), "selection", net * "_thresholds.csv")

"Stability class of each ParamID = number of distinct discrete states among its steady states."
function stability_classes(rows, θ)
    states = Dict{Int,Set{Int}}()
    for (pid, _, logx) in rows
        push!(get!(states, pid, Set{Int}()), state_code(2.0 .^ logx, θ))
    end
    return Dict(pid => get(CLASS_NAMES, length(s), "multistable_$(length(s))") for (pid, s) in states)
end

function select_parameters(cfg::Config, netname::AbstractString)
    verify_inputs(cfg, netname)
    net = load_network(cfg, netname)
    rows = read_solutions(cfg, net)
    θ = thresholds(rows, nnodes(net))
    cls = stability_classes(rows, θ)
    params = read_parameters(cfg, net)
    keys(cls) == keys(params) || error("$netname: solution and parameter files list different ParamIDs")
    ntotal = length(cls)

    out = DataFrame(ParamID=Int[], StabilityClass=String[], ClassSize=Int[],
                    ClassPrevalence=Float64[], NSelected=Int[], SelectionProb=Float64[])
    for c in cfg.classes
        ids = sort([pid for (pid, k) in cls if k == c])
        isempty(ids) && continue
        n = min(cfg.n_per_class, length(ids))
        rng = rng_for(cfg.seed_root, netname, "selection", c)
        chosen = sort(sample(rng, ids, n; replace=false))
        for pid in chosen
            push!(out, (pid, c, length(ids), length(ids) / ntotal, n, n / length(ids)))
        end
    end
    sort!(out, :ParamID)

    mkpath(dirname(selection_path(cfg, netname)))
    CSV.write(selection_path(cfg, netname), out)
    CSV.write(thresholds_path(cfg, netname), DataFrame(Node=net.nodes, Threshold=θ))
    open(selection_path(cfg, netname) * ".meta.toml", "w") do io
        TOML.print(io, Dict("config_hash" => cfg.hash, "network" => netname,
                            "n_param_sets_total" => ntotal, provenance()...))
    end
    return out
end

"Load a selection, refusing one made with a different config."
function load_selection(cfg::Config, netname::AbstractString)
    path = selection_path(cfg, netname)
    isfile(path) || error("no selection for $netname at $path — run scripts/select.jl first")
    meta = TOML.parsefile(path * ".meta.toml")
    meta["config_hash"] == cfg.hash ||
        error("selection for $netname was made with config $(meta["config_hash"]), current is $(cfg.hash)")
    # Plain String columns: newer CSV.jl returns its own string types, which would then
    # leak into the raw files and need that exact package to be read back.
    sel = CSV.read(path, DataFrame; types=Dict(:StabilityClass => String))
    sel.StabilityClass = String.(sel.StabilityClass)
    th = CSV.read(thresholds_path(cfg, netname), DataFrame; types=Dict(:Node => String))
    return sel, Vector{Float64}(th.Threshold), String.(th.Node)
end
