# Minimal --key value parser shared by the scripts. Unknown flags are errors.
function parse_cli(args, allowed::Dict{String,Any})
    out = copy(allowed)
    i = 1
    while i <= length(args)
        a = args[i]
        startswith(a, "--") || error("unexpected argument '$a'")
        k = a[3:end]
        haskey(allowed, k) || error("unknown flag --$k (allowed: $(join(sort(collect(keys(allowed))), ", ")))")
        if allowed[k] isa Bool
            out[k] = true; i += 1
        else
            i < length(args) || error("--$k needs a value")
            out[k] = args[i+1]; i += 2
        end
    end
    return out
end

const PKG = normpath(joinpath(@__DIR__, ".."))
default_config() = joinpath(PKG, "config", "production.toml")

cfg_from(o) = RegNoise.load_config(o["config"];
    data_dir=o["data-dir"] === nothing ? nothing : o["data-dir"],
    results_root=o["results-root"] === nothing ? nothing : o["results-root"])

nets_from(o, cfg) = o["network"] === nothing ? cfg.networks : split(o["network"], ",")
modes_from(o, cfg) = o["mode"] === nothing ? cfg.noise_modes : split(o["mode"], ",")
