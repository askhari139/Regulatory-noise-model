using JLD2

"""
    run_parallel_shell(param_ids, noise_mode, param_type, num_sims, noise_levels, stability_class; 
                      max_parallel=4, log_dir="shell_logs")

Run multiple parameter IDs in parallel using shell background jobs
"""

"""
    run_parallel_shell(param_ids, noise_mode, param_type, num_sims, noise_levels, stability_class; 
                      max_parallel=4, log_dir="shell_logs")
 
Run multiple parameter IDs in parallel using shell background jobs
"""
 
function run_parallel_shell(param_ids::Vector{Int},
                           noise_mode::String,
                           network::String,
                           param_type::Symbol,
                           num_sims::Int,
                           noise_levels::Vector{Float64},
                           stability_class::Vector{Bool};
                           max_parallel::Int=4,
                           log_dir::String="shell_logs",
                           work_dir::String=joinpath(noise_mode, network),
                           deterministic::Bool=false,
                           resample::Bool=false,
                           det_iters::Int=10,
                           source_dir::String=sourceDir)
 
    mkpath(log_dir)
    println("Precompiling packages (one-time overhead)...")
    run(`julia --project=. -e 'using Pkg; Pkg.precompile()'`)
    println("✓ Precompilation complete\n")
 
    noise_levels_str = join(noise_levels, ",")
 
    println("="^70)
    println("SHELL PARALLEL MODE")
    println("="^70)
    println("Running $(length(param_ids)) parameters")
    println("Max parallel jobs: $max_parallel")
    println("Deterministic: $deterministic")
    println("="^70)
 
    script_name = deterministic ? "3_analyze_single_parameter_det.jl" :
                                  "3_analyze_single_parameter.jl"
    script = joinpath(source_dir, "scripts", script_name)
    isfile(script) || error("Script not found: $script")
 
    n_batches = ceil(Int, length(param_ids) / max_parallel)
 
    for batch_num in 1:n_batches
        batch_start = (batch_num - 1) * max_parallel + 1
        batch_end   = min(batch_num * max_parallel, length(param_ids))
        batch_params = param_ids[batch_start:batch_end]
 
        println("\nBatch $batch_num/$n_batches: parameters $batch_params")
 
        processes = Tuple{Int, Base.Process}[]
        for param_id in batch_params
            logfile = joinpath(log_dir, "param_$(param_id).log")
 
            cmd = `julia --project=. --threads=1
                   $script
                   param-id=$(string(param_id))
                   noise-mode=$(string(noise_mode))
                   network=$network
                   work-dir=$work_dir
                   noise-levels=$noise_levels_str
                   resample=$(string(resample))
                   det-iters=$det_iters`
 
            proc = run(pipeline(cmd, stdout=logfile, stderr=logfile), wait=false)
            push!(processes, (param_id, proc))
            try
                println("  Started parameter $param_id (PID: $(getpid(proc)))")
            catch
                println("  Started parameter $param_id")
            end
        end
 
        println("\n  Waiting for batch to complete...")
        for (param_id, proc) in processes
            wait(proc)
            if success(proc)
                println("  ✓ Parameter $param_id completed")
            else
                @warn "  ✗ Parameter $param_id failed — check $log_dir/param_$(param_id).log"
            end
        end
    end
end


function submit_slurm_job_array(param_ids, noise_mode, network, param_type, num_sims,
                               noise_levels, stability_class, log_dir;
                               sourceDir = ".",
                               deterministic::Bool = false, #use_uniform::Bool = false,
                               resample::Bool = false,
                               cpus=1, mem="8G", time="2:00:00", partition="ctbp", det_iters::Int=10)
    
    mkpath(log_dir)
    
    println("Precompiling packages (one-time)...")
    run(`julia --project=. -e 'using Pkg; Pkg.precompile()'`)
    println("✓ Precompilation complete\n")
    
    param_file = joinpath(log_dir, "param_ids.txt")
    open(param_file, "w") do f
        for param_id in param_ids
            println(f, param_id)
        end
    end
    
    noise_levels_str = join(noise_levels, ",")
    stability_str = join([x ? "T" : "F" for x in stability_class], "")
    SCRIPT = deterministic ? "3_analyze_single_parameter_det.jl" : "3_analyze_single_parameter.jl"
    
    # Detect cluster max array size
    max_array_size = 1000
    try
        config_out = read(`scontrol show config`, String)
        m = match(r"MaxArraySize\s*=\s*(\d+)", config_out)
        if m !== nothing
            max_array_size = parse(Int, m.captures[1]) - 1
        end
    catch
        @warn "Could not query MaxArraySize, defaulting to $max_array_size"
    end
    
    n_total  = length(param_ids)
    n_chunks = ceil(Int, n_total / max_array_size)
    
    println("Submitting $n_total tasks in $n_chunks array(s) (max size: $max_array_size)")
    
    # source_dir = get(ENV, "RACIPE_SOURCE", pwd())
    source_dir = sourceDir
    work_dir = joinpath(noise_mode, network)
    
    for chunk in 1:n_chunks
        i_start = (chunk - 1) * max_array_size + 1
        i_end   = min(chunk * max_array_size, n_total)
        chunk_size = i_end - i_start + 1
        
        job_script = """
#!/bin/bash
#SBATCH --job-name=stoch_$(noise_mode)_$(chunk)
#SBATCH --output=$log_dir/param_%a.out
#SBATCH --error=$log_dir/param_%a.err
#SBATCH --array=1-$chunk_size
#SBATCH --cpus-per-task=$cpus
#SBATCH --partition=$partition
#SBATCH --mem=$mem
#SBATCH --time=$time

set -euo pipefail

cd $(pwd())

echo "=== STARTUP ==="
echo "Host:    \$(hostname)"
echo "PWD:     \$(pwd)"
echo "Chunk:   $chunk / $n_chunks"
echo "Task ID: \${SLURM_ARRAY_TASK_ID}"

# Offset task ID to global param index
GLOBAL_IDX=\$(( \${SLURM_ARRAY_TASK_ID} + $(i_start - 1) ))
PARAM_ID=\$(sed -n "\${GLOBAL_IDX}p" $param_file)
echo "Global:  \${GLOBAL_IDX}"
echo "Param:   \${PARAM_ID}"
echo "==============="

export RACIPE_SOURCE=$source_dir

echo "Julia Starts."
echo "HOME: \$HOME"
ls ~/.julia/compiled/ > /dev/null 2>&1 && echo "depot OK" || echo "depot FAIL"
stdbuf -oL -eL julia --project=. --threads=$cpus \\
      $source_dir/scripts/$SCRIPT \\
      noise-mode=$noise_mode \\
      network=$network \\
      param-type=$(string(param_type)) \\
      num-sims=$num_sims \\
      noise-levels=$noise_levels_str \\
      param-id=\${PARAM_ID} \\
      stability-class=$stability_str \\
      work-dir=$work_dir \\
      resample=$(string(resample)) \\
      det-iters=$det_iters
echo "Julia ends."
"""
        
        chunk_file = joinpath(log_dir, "submit_array_$(chunk).sh")
        open(chunk_file, "w") do f
            write(f, job_script)
        end
        
        run(`sbatch $chunk_file`)
        println("  Submitted chunk $chunk/$n_chunks (tasks $i_start – $i_end)")
    end
    
    println("\n✓ Submitted $n_total tasks across $n_chunks job array(s)")
    println("Monitor with: squeue -u \$USER")
end

function run_parallel_cluster(param_ids::Vector{Int}, 
                           noise_mode::String,
                           param_type::Symbol,
                           num_sims::Int,
                           noise_levels::Vector{Float64},
                           stability_class::Vector{Bool};
                           log_dir::String="cluster_logs")
    mkpath(log_dir)
    noise_levels_str = join(noise_levels, ",")
    stability_str = join([x ? "T" : "F" for x in stability_class], "")
    script = joinpath(scriptsDir, 3_analyze_single_parameter.jl)
    for param_id in param_ids
        cmd = ```
                julia --project=. --threads=1 
                $script 
                param-id=$(string(param_id))
                ```
        job_script = """
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=param_$(param_id)
#SBATCH --output=cluster_logs/param_$(param_id).out
#SBATCH --error=cluster_logs/param_$(param_id).err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00

cd $(pwd())
$(cmd)
"""
            
        # Write job script
        job_file = "cluster_logs/submit_$(param_id).sh"
        open(job_file, "w") do f
            write(f, job_script)
        end
        
        # Submit job
        run(`sbatch $job_file`)
        println("  Submitted parameter $param_id ($param_type)")
        
        # Small delay to avoid overwhelming scheduler
        sleep(0.1)
    end
    println("\n✓ Submitted $(length(param_ids)) jobs")
    println("Monitor with: squeue -u \$USER")
    println("Check logs in: cluster_logs/")
    println("\nTo collect results after jobs finish, run:")
    println("  julia --project=. scripts/collect_results.jl $net")
end

