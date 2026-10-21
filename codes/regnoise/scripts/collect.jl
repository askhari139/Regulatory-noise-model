# Validate, then merge raw unit files into gzipped long-format CSVs (SPEC §11).
# julia --project scripts/collect.jl [--network TS] [--mode Additive]
using RegNoise, CSV, DataFrames, CodecZlib, Dates
include(joinpath(@__DIR__, "cli.jl"))
include(joinpath(@__DIR__, "validation.jl"))
o = parse_cli(ARGS, Dict{String,Any}("config" => default_config(), "network" => nothing,
    "mode" => nothing, "data-dir" => nothing, "results-root" => nothing))
cfg = cfg_from(o)
for netname in nets_from(o, cfg), mode in modes_from(o, cfg)
    validate_network_mode(cfg, netname, mode) ||
        error("$netname $mode failed validation; not collecting (see validation/$(netname)_$(mode).txt)")
    sel, _, _ = load_selection(cfg, netname)
    outdir = joinpath(run_dir(cfg), "collected"); mkpath(outdir)
    streams = Dict{String,Any}()
    started = Set{String}()
    try
        for pid in sel.ParamID
            meta, T = RegNoise.read_unit(raw_path(cfg, netname, mode, pid))
            for (tname, df) in T
                nrow(df) == 0 && continue
                io = get!(streams, tname) do
                    GzipCompressorStream(open(joinpath(outdir, "$(netname)_$(mode)_$(tname).csv.gz"), "w"))
                end
                df = insertcols!(copy(df), 1, :RunID => cfg.run_id, :ConfigHash => cfg.hash,
                                 :Network => netname, :NoiseMode => mode, :ParamID => pid,
                                 :StabilityClass => meta["stability_class"])
                CSV.write(io, df; append=tname in started, writeheader=!(tname in started))
                push!(started, tname)
            end
        end
    finally
        foreach(close, values(streams))
    end
    println("collected $netname $mode -> $outdir")
end
