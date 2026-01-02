#!/usr/bin/env julia
"""
Estimate build times for Pluto notebooks to help plan CI parallelization.

Usage:
    julia services/estimate_build_time.jl [--sample N] [--timeout SECONDS]
    
Options:
    --sample N        Number of notebooks to sample (default: 3)
    --timeout SECONDS Maximum time per notebook in seconds (default: 300)
    
This script samples a few notebooks, measures their build time, and 
extrapolates to estimate total build time for sequential vs parallel builds.
"""

using Pkg

# Parse command line arguments
function parse_args()
    sample_count = 3
    timeout = 300
    
    args = ARGS
    i = 1
    while i <= length(args)
        if args[i] == "--sample" && i < length(args)
            sample_count = parse(Int, args[i+1])
            i += 2
        elseif args[i] == "--timeout" && i < length(args)
            timeout = parse(Int, args[i+1])
            i += 2
        else
            i += 1
        end
    end
    
    return (; sample_count, timeout)
end

function find_pluto_notebooks(dir::String)
    notebooks = String[]
    for file in readdir(dir; join=true)
        if endswith(file, ".jl")
            try
                first_line = readline(file)
                if first_line == "### A Pluto.jl notebook ###"
                    push!(notebooks, file)
                end
            catch
            end
        end
    end
    return notebooks
end

function main()
    opts = parse_args()
    
    println("="^60)
    println("Pluto Notebook Build Time Estimator")
    println("="^60)
    
    # Find notebooks
    notebooks_dir = joinpath(@__DIR__, "..", "notebooks")
    all_notebooks = find_pluto_notebooks(notebooks_dir)
    
    println("\nFound $(length(all_notebooks)) Pluto notebooks")
    
    if isempty(all_notebooks)
        println("No notebooks found!")
        return
    end
    
    # Sample notebooks for timing
    sample_size = min(opts.sample_count, length(all_notebooks))
    
    # Choose diverse samples (first, middle, last to get variety)
    indices = unique([1, length(all_notebooks) ÷ 2, length(all_notebooks)])
    if length(indices) < sample_size
        # Add random samples
        remaining = setdiff(1:length(all_notebooks), indices)
        additional = rand(remaining, sample_size - length(indices))
        indices = sort(unique([indices..., additional...]))
    end
    indices = indices[1:min(sample_size, length(indices))]
    
    sample_notebooks = all_notebooks[indices]
    
    println("\nSampling $(length(sample_notebooks)) notebooks for timing:")
    for nb in sample_notebooks
        println("  - $(basename(nb))")
    end
    
    # Install PlutoSliderServer if needed
    println("\nSetting up PlutoSliderServer...")
    tmp_env = mktempdir()
    Pkg.activate(tmp_env)
    Pkg.add(Pkg.PackageSpec(name="PlutoSliderServer", version="1"); io=devnull)
    
    using PlutoSliderServer: export_notebook
    
    # Time each notebook
    println("\n" * "="^60)
    println("Timing notebooks (this may take several minutes)...")
    println("="^60)
    
    timings = Dict{String, Float64}()
    output_dir = mktempdir()
    cache_dir = mktempdir()
    
    for notebook in sample_notebooks
        name = basename(notebook)
        println("\nBuilding: $name")
        
        try
            elapsed = @elapsed begin
                export_notebook(
                    notebook;
                    Export_cache_dir=cache_dir,
                    Export_baked_notebookfile=false,
                    Export_baked_state=false,
                    Export_output_dir=output_dir,
                )
            end
            timings[name] = elapsed
            println("  ✓ Completed in $(round(elapsed, digits=1))s")
        catch e
            println("  ✗ Failed: $e")
            timings[name] = NaN
        end
    end
    
    # Calculate statistics
    valid_timings = filter(!isnan, collect(values(timings)))
    
    if isempty(valid_timings)
        println("\nNo successful builds to analyze!")
        return
    end
    
    avg_time = sum(valid_timings) / length(valid_timings)
    max_time = maximum(valid_timings)
    min_time = minimum(valid_timings)
    total_notebooks = length(all_notebooks)
    
    println("\n" * "="^60)
    println("TIMING ANALYSIS")
    println("="^60)
    
    println("\nSample Results:")
    println("  - Average time per notebook: $(round(avg_time, digits=1))s")
    println("  - Min time: $(round(min_time, digits=1))s")
    println("  - Max time: $(round(max_time, digits=1))s")
    
    println("\nExtrapolation to all $total_notebooks notebooks:")
    
    # Sequential estimate
    seq_estimate = avg_time * total_notebooks
    println("\n  SEQUENTIAL BUILD (current approach):")
    println("    Estimated time: $(round(seq_estimate / 60, digits=1)) minutes")
    println("    ($(round(seq_estimate, digits=0))s)")
    
    # Parallel estimates for different chunk sizes
    println("\n  PARALLEL BUILD OPTIONS:")
    
    for chunk_size in [3, 5, 10]
        num_chunks = ceil(Int, total_notebooks / chunk_size)
        # Parallel time ≈ max chunk time + overhead
        # Assume max chunk takes avg_time * chunk_size with some variance
        parallel_estimate = avg_time * chunk_size * 1.2  # 20% overhead for variance
        speedup = seq_estimate / parallel_estimate
        
        println("\n    $(num_chunks) parallel jobs ($chunk_size notebooks each):")
        println("      Estimated time: $(round(parallel_estimate / 60, digits=1)) minutes")
        println("      Speedup: ~$(round(speedup, digits=1))x")
    end
    
    println("\n" * "="^60)
    println("RECOMMENDATIONS")
    println("="^60)
    
    if seq_estimate > 3600  # More than 1 hour
        println("\n⚠️  Sequential build would take over 1 hour!")
        println("   STRONGLY recommend parallel builds with 5-10 chunks")
    elseif seq_estimate > 1800  # More than 30 min
        println("\n⚡ Sequential build would take 30-60 minutes")
        println("   Recommend parallel builds with 5 chunks")
    else
        println("\n✓ Build times are manageable")
        println("   Parallel builds would still provide ~5x speedup")
    end
    
    println()
end

main()
