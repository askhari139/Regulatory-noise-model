"""
Module to run RACIPE and manage output files
"""

module RACIPErunner

export run_racipe, run_racipe_batch, clean_racipe_output

"""
    run_racipe(topo_file; num_paras=10000, threads=4, keep_only_essentials=true)

Run RACIPE on a topology file and optionally clean up extra files

# Arguments
- `topo_file`: Path to .topo file (e.g., "TS.topo" or "data/TS.topo")
- `num_paras`: Number of parameter sets to generate (default: 10000)
- `threads`: Number of threads to use (default: 4)
- `keep_only_essentials`: If true, keep only .prs, _solution.dat, _parameters.dat (default: true)

# Returns
- Tuple of (prs_file, solution_file, parameters_file) paths

# Example
```julia
prs, sol, params = run_racipe("data/TS.topo", num_paras=5000, threads=8)
```
"""
function run_racipe(topo_file::String; 
                    num_paras::Int=10000, 
                    threads::Int=4,
                    keep_only_essentials::Bool=true,
                    force::Bool=false)
    
    # Check if topology file exists
    if !isfile(topo_file)
        error("Topology file not found: $topo_file")
    end
    
    # Check if racipemt is available
    if !check_racipe_installed()
        error("racipemt not found in PATH. Please install RACIPE.")
    end
    
    # Get base name without extension
    base_name = replace(basename(topo_file), r"\.topo$" => "")
    topo_dir = dirname(topo_file)
    if isempty(topo_dir)
        topo_dir = "."
    end
    runRacipe = false
    # Expected output files
    prs_file = joinpath(topo_dir, "$base_name.prs")
    solution_file = joinpath(topo_dir, "$(base_name)_solution.dat")
    parameters_file = joinpath(topo_dir, "$(base_name)_parameters.dat")

    # Check if essential files are already present
    essential_files = [prs_file, solution_file, parameters_file]
    missing_files = [f for f in essential_files if !isfile(f)]
    
    if !isempty(missing_files) || force
        runRacipe = true
    else
        n_lines = countlines(parameters_file)
        if n_lines < num_paras
            runRacipe = true
        end
    end
    if !runRacipe
        println("A run of racipe already exists with the number of parameters requested. If a fresh run is required, set force = true.")
        return (prs_file, solution_file, parameters_file)
    end
    if runRacipe
        println("="^60)
        println("Running RACIPE")
        println("="^60)
        println("Topology file: $topo_file")
        println("Base name: $base_name")
        println("Parameters: $num_paras")
        println("Threads: $threads")
        println("Output directory: $topo_dir")
        println()
        
        mkpath("TempSim")
        cp(topo_file, "TempSim/"*topo_file)
        cd("TempSim")
        # Build command
        cmd = `racipemt $topo_file -threads $threads -num_paras $num_paras`
        
        # Run RACIPE
        println("Running command: $cmd")
        println("This may take a few minutes...")
        println()
        
        try
            run(cmd)
            println("\n✓ RACIPE completed successfully!")
        catch e
            error("RACIPE execution failed: $e")
        end
        
        missing_files = [f for f in essential_files if !isfile(f)]
        
        if !isempty(missing_files)
            error("Expected output files not found: $(join(missing_files, ", "))")
        end
        
        println("\nEssential output files created:")
        println("  ✓ $prs_file")
        println("  ✓ $solution_file")
        println("  ✓ $parameters_file")
    end
    
    # Clean up extra files if requested
    if keep_only_essentials
        println("\nCleaning up extra output files...")
        deleted_count = clean_racipe_output(topo_dir, base_name)
        if deleted_count > 0
            println("  ✓ Deleted $deleted_count extra file(s)")
        else
            println("  ✓ No extra files to delete")
        end
    end
    
    println("\n" * "="^60)
    println("RACIPE run complete!")
    println("="^60)
    cp(prs_file, "../$prs_file")
    cp(solution_file, "../$solution_file")
    cp(parameters_file, "../$parameters_file")
    cd("..")
    rm("TempSim", recursive=true, force=true)
    return (prs_file, solution_file, parameters_file)
end

"""
    check_racipe_installed()

Check if racipemt is available in PATH
"""
function check_racipe_installed()
    try
        run(pipeline(`which racipemt`, stdout=devnull, stderr=devnull))
        return true
    catch
        try
            # Try running it directly
            run(pipeline(`racipemt`, stdout=devnull, stderr=devnull))
            return true
        catch
            return false
        end
    end
end

"""
    clean_racipe_output(directory, base_name)

Delete all RACIPE output files except .prs, _solution.dat, and _parameters.dat

Returns: Number of files deleted
"""
function clean_racipe_output(directory::String, base_name::String)
    # Files to keep
    keep_files = Set([
        "$base_name.prs",
        "$(base_name)_solution.dat",
        "$(base_name)_parameters.dat",
        "$(base_name).topo"
    ])
    
    deleted_count = 0
    
    # List all files in directory
    for file in readdir(directory)
        # Skip if it's a file we want to keep
        if file in keep_files
            continue
        else
            file_path = joinpath(directory, file)
            try
                rm(file_path)
                println("    Deleted: $file")
                deleted_count += 1
            catch e
                @warn "Could not delete $file: $e"
            end
        end
    end
    
    return deleted_count
end

"""
    run_racipe_batch(topo_files; num_paras=10000, threads=4, keep_only_essentials=true)

Run RACIPE on multiple topology files

# Arguments
- `topo_files`: Vector of paths to .topo files
- Other arguments same as run_racipe

# Returns
- Dictionary mapping base_name => (prs_file, solution_file, parameters_file)

# Example
```julia
results = run_racipe_batch(["TS.topo", "EMT.topo"], num_paras=5000)
```
"""
function run_racipe_batch(topo_files::Vector{String}; 
                         num_paras::Int=10000, 
                         threads::Int=4,
                         keep_only_essentials::Bool=true)
    
    results = Dict{String, Tuple{String, String, String}}()
    
    println("\n" * "="^60)
    println("BATCH RACIPE EXECUTION")
    println("="^60)
    println("Number of networks: $(length(topo_files))")
    println("Parameters per network: $num_paras")
    println("Threads: $threads")
    println("="^60)
    println()
    
    for (i, topo_file) in enumerate(topo_files)
        println("\n[$i/$(length(topo_files))] Processing: $topo_file")
        println("-"^60)
        
        try
            base_name = replace(basename(topo_file), r"\.topo$" => "")
            files = run_racipe(topo_file, 
                             num_paras=num_paras, 
                             threads=threads,
                             keep_only_essentials=keep_only_essentials)
            results[base_name] = files
            println("✓ Success: $base_name")
        catch e
            @error "Failed to process $topo_file: $e"
        end
    end
    
    println("\n" * "="^60)
    println("BATCH COMPLETE")
    println("="^60)
    println("Successfully processed: $(length(results))/$(length(topo_files)) networks")
    
    return results
end

end # module
