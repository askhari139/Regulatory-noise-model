"""
    parse_named_args(args=ARGS; defaults=Dict())

Parse command line arguments in key=value format with defaults

# # Example
# ```julia
# args = parse_named_args(defaults=Dict("mode" => "serial", "n" => "10"))
# mode = get_arg(args, "mode", "serial")
# n = get_arg(args, "n", 10, type=Int)
# ```
# """
function parse_named_args(args=ARGS; defaults=Dict())
    parsed = copy(defaults)
    
    for arg in args
        if occursin("=", arg)
            # Handle --key=value or key=value
            arg = replace(arg, r"^--" => "")
            key, val = split(arg, "=", limit=2)
            parsed[key] = val
        end
    end
    
    return parsed
end

"""
    split_nonempty(val, sep)

Split on `sep` and drop blank fields (after stripping whitespace) -- plain
`split` keeps empty fields, so a trailing/leading/doubled separator (e.g.
"TTSA," or "TS,,TT") silently produces an empty-string element downstream
(e.g. an empty "network" that then fails on a ".topo" file lookup) instead
of erroring loudly or being ignored.
"""
function split_nonempty(val, sep::String)
    return [string(strip(x)) for x in split(string(val), sep) if !isempty(strip(x))]
end

"""
    get_arg(args, key, default; type=String)

Get argument value with type conversion
"""
function get_arg(args::Dict, key::String, default; type::Type=String, sep::String = ",")
    
    val = get(args, key, default)
    # println(val)
    val = string(val)
    if type == Int
        return parse(Int, val)
    elseif type == Float64
        return parse(Float64, val)
    elseif type == Symbol
        return Symbol(val)
    elseif type == Bool
        return lowercase(val) in ["true", "1", "yes"]
    elseif type == Vector{Float64}
        return parse.(Float64, split_nonempty(val, sep))
    elseif type == Vector{Int}
        return parse.(Int, split_nonempty(val, sep))
    elseif type == Vector{String}
        return split_nonempty(val, sep)
    elseif type == Vector{Bool}
        if typeof(val) == Vector{Bool}
            return val
        elseif occursin(sep, string(val))
            return [lowercase(string(x)) in ["true", "1", "yes", "t"] for x in split_nonempty(string(val), sep)]
        else
            # String like "TTFF"
            return [c in ['T', 't', '1'] for c in string(val)]
        end
    else
        return string(val)
    end
end

"""
    print_config(config::Dict)

Pretty-print configuration
"""
function print_config(config::Dict)
    println("="^70)
    println("CONFIGURATION")
    println("="^70)
    for (key, val) in sort(collect(config))
        println("  $key = $val")
    end
    println("="^70)
end