# Shared by validate.jl and collect.jl.
function validate_network_mode(cfg, netname, mode)
    net = load_network(cfg, netname)
    sel, _, _ = load_selection(cfg, netname)
    lo, hi = lambda_bounds(cfg, net)
    report = String[]
    missing_ids = Int[]
    for pid in sel.ParamID
        path = raw_path(cfg, netname, mode, pid)
        if !isfile(path)
            push!(missing_ids, pid); continue
        end
        for p in validate_unit(cfg, path, mode, lo, hi)
            push!(report, "ParamID $pid: $p")
        end
    end
    isempty(missing_ids) || push!(report, "$(length(missing_ids)) of $(nrow(sel)) ParamIDs missing: $(join(first(missing_ids, 20), ","))$(length(missing_ids) > 20 ? ",..." : "")")
    out = joinpath(run_dir(cfg), "validation", "$(netname)_$(mode).txt")
    mkpath(dirname(out))
    status = isempty(report) ? "PASS" : "FAIL"
    open(out, "w") do io
        println(io, "$status  $netname $mode  run $(cfg.run_id)  config $(cfg.hash)  $(now())")
        println(io, "$(nrow(sel)) selected ParamIDs; expected cells per unit: $(length(expected_cells(cfg, mode)))")
        foreach(l -> println(io, l), report)
    end
    println("$status  $netname $mode  ($(length(report)) problems) -> $out")
    return isempty(report)
end
