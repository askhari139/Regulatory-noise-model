# Select parameter sets for each network (SPEC §4). Run once per run_id.
# julia --project scripts/select.jl [--config FILE] [--network TS,TT]
using RegNoise
include(joinpath(@__DIR__, "cli.jl"))
o = parse_cli(ARGS, Dict{String,Any}("config" => default_config(), "network" => nothing,
    "data-dir" => nothing, "results-root" => nothing, "allow-dirty" => false,
    "skip-existing" => false))
cfg = cfg_from(o)
check_clean_tree(o["allow-dirty"])
for net in nets_from(o, cfg)
    if o["skip-existing"] && isfile(RegNoise.selection_path(cfg, net))
        load_selection(cfg, net)   # errors if it was made with another config
        println("$net: selection exists, keeping it"); continue
    end
    sel = select_parameters(cfg, net)
    counts = [(c, count(==(c), sel.StabilityClass)) for c in unique(sel.StabilityClass)]
    println("$net: ", join(["$c=$n" for (c, n) in counts], ", "), "  (config $(cfg.hash[1:12]))")
end
