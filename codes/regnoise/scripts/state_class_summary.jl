# Sum occupancy over state classes (by number of high nodes) for every parameter set and
# condition cell, from collected occupancy files. Output feeds the R figure scripts.
#   julia --project scripts/state_class_summary.jl --collected DIR --out FILE.csv
# Classes: all-low (0 high), single-high (1), double-high (2), triple-high (3), all-high (N).
using CSV, DataFrames, CodecZlib
include(joinpath(@__DIR__, "cli.jl"))
o = parse_cli(ARGS, Dict{String,Any}("collected" => nothing, "out" => nothing))
(o["collected"] === nothing || o["out"] === nothing) && error("--collected and --out are required")

function class_name(label)
    k = count(==('H'), label); n = length(label)
    k == n ? "all-high" : k == 0 ? "all-low" : ("single-high", "double-high", "triple-high")[k]
end

out = DataFrame()
for f in sort(readdir(o["collected"]))
    endswith(f, "_occupancy.csv.gz") || continue
    # stringtype=String: CSV.jl otherwise picks fixed-width string types per file (String3 for
    # "TS"), which cannot hold longer values from the next file ("TSSA") when appending
    occ = CSV.read(GzipDecompressorStream(open(joinpath(o["collected"], f))), DataFrame;
                   stringtype=String, types=Dict(:State => String, :StabilityClass => String))
    occ = filter(r -> r.Block == 0 && r.ThresholdMult == 1.0, occ)
    occ.StateClass = class_name.(occ.State)
    g = combine(groupby(occ, [:Network, :NoiseMode, :ParamID, :StabilityClass, :Condition,
                              :Sigma, :DT, :THold, :StateClass]), :Mean => sum => :Occupancy)
    append!(out, g)
    println("$f: $(nrow(g)) rows")
end
CSV.write(o["out"], out)
println("wrote $(o["out"])")
