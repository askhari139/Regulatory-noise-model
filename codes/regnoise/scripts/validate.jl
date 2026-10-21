# Validate raw outputs against the config (SPEC §11). Exit code 1 on any problem.
# julia --project scripts/validate.jl [--network TS] [--mode Additive]
using RegNoise, DataFrames, Dates
include(joinpath(@__DIR__, "cli.jl"))
include(joinpath(@__DIR__, "validation.jl"))
o = parse_cli(ARGS, Dict{String,Any}("config" => default_config(), "network" => nothing,
    "mode" => nothing, "data-dir" => nothing, "results-root" => nothing))
cfg = cfg_from(o)
ok = true
for net in nets_from(o, cfg), mode in modes_from(o, cfg)
    global ok &= validate_network_mode(cfg, net, mode)
end
ok || exit(1)
