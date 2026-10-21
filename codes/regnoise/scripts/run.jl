# Run all conditions for a chunk of selected ParamIDs (SPEC §6, §13).
# julia --project scripts/run.jl --network TS --mode Additive [--chunk 3/40] [--config FILE]
using RegNoise
include(joinpath(@__DIR__, "cli.jl"))
o = parse_cli(ARGS, Dict{String,Any}("config" => default_config(), "network" => nothing,
    "mode" => nothing, "chunk" => "1/1", "params" => nothing,
    "data-dir" => nothing, "results-root" => nothing, "allow-dirty" => false))
o["network"] === nothing && error("--network is required")
o["mode"] === nothing && error("--mode is required")
cfg = cfg_from(o)
check_clean_tree(o["allow-dirty"])
netname, mode = o["network"], o["mode"]
netname in cfg.networks || error("network $netname is not in the config")
mode in cfg.noise_modes || error("mode $mode is not in the config")

verify_inputs(cfg, netname)
net = load_network(cfg, netname)
sel, θ, nodes = load_selection(cfg, netname)
nodes == net.nodes || error("selection node order $nodes differs from network $(net.nodes)")
params = read_parameters(cfg, net)

k, K = parse.(Int, split(o["chunk"], "/"))
1 <= k <= K || error("bad --chunk $(o["chunk"])")
ids = o["params"] === nothing ? sel.ParamID[k:K:end] : parse.(Int, split(o["params"], ","))
class = Dict(zip(sel.ParamID, sel.StabilityClass))
for pid in ids
    haskey(class, pid) || error("ParamID $pid is not in the selection for $netname")
end

println("run $(cfg.run_id) | $netname $mode | chunk $k/$K | $(length(ids)) param sets | config $(cfg.hash[1:12])")
failed = Int[]
for pid in ids
    if unit_status(cfg, netname, mode, pid) == :done
        println("  $pid: done, skipping"); continue
    end
    try
        t = @elapsed rows = run_unit(cfg, net, mode, pid, params[pid], θ)
        write_unit(cfg, net, mode, pid, class[pid], θ, rows, t)
        println("  $pid ($(class[pid])): $(round(t, digits=1)) s")
    catch e
        push!(failed, pid)
        println("  $pid: FAILED")
        showerror(stdout, e, catch_backtrace()); println()
    end
    flush(stdout)
end
isempty(failed) || (println("FAILED ParamIDs: ", join(failed, ",")); exit(1))
