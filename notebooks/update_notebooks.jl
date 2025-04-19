using Pluto: Pluto
using Pkg: Pkg
import StatsBase: countmap
using Logging

"""
Update all Pluto notebooks in the current directory and check version compatibility.
This script:
1. Updates all notebooks starting with "N-0"
2. Upgrades manifests to v2.0.0 if needed
3. Checks for conflicting major versions of packages across notebooks
"""

# Configure logging with more detail
logger = ConsoleLogger(stderr, Logging.Debug)
global_logger(logger)

# Track failures for summary
failed_notebooks = Dict{String, Any}()
successful_updates = String[]

# Find all notebooks matching pattern
notebooks = readdir(@__DIR__)
filter!(notebooks) do name
    contains(name, r"N-0")
end

@info "Found $(length(notebooks)) notebooks to process"

# Update each notebook's environment
for name in notebooks
    @info "Processing notebook: $name"
    try
        notebook_path = joinpath(@__DIR__, name)
        @debug "Full path: $notebook_path"
        
        if !isfile(notebook_path)
            throw(ErrorException("Notebook file not found"))
        end
        
        Pluto.activate_notebook_environment(notebook_path)
        
        # Check and upgrade manifest if needed
        ctx = Pkg.API.Context()
        if ctx.env.manifest.manifest_format != v"2.0.0"
            @info "Upgrading manifest for $name"
            Pkg.upgrade_manifest()
            Pkg.resolve()
        end
        
        # Update packages
        Pkg.update()
        push!(successful_updates, name)
        
    catch e
        @error "Failed to process notebook $name" exception=(e, catch_backtrace())
        failed_notebooks[name] = (error=e, trace=catch_backtrace())
    end
end

# Restore original environment
Pkg.activate(joinpath(@__DIR__, ".."))

"""
Extract compatibility requirements from a Pluto notebook file.
Returns a tuple of (success::Bool, compat::Vector{String}, error::Union{Nothing,Exception})
"""
function getcompat(file)
    compat = String[]
    start_fill = false
    
    try
        if !isfile(file)
            return (false, String[], ErrorException("File not found: $file"))
        end
        
        for line in eachline(file)
            if line == "[compat]"
                start_fill = true
                continue
            end
            if start_fill
                if line == "\"\"\""
                    break
                end
                push!(compat, line)
            end
        end
        return (true, compat, nothing)
    catch e
        return (false, String[], e)
    end
end

"""
Parse a version specification line from the [compat] section.
Returns a tuple of (package_name, version) or nothing if parsing fails.
Example input: "Package = \"1.0.0\""
"""
function parse_version_line(line)
    try
        # Remove any whitespace and quotes
        line = strip(line)
        # Skip empty lines or comments
        if isempty(line) || startswith(line, "#")
            return nothing
        end
        
        # Split on equals sign and clean up parts
        parts = split(line, "=", limit=2)
        if length(parts) != 2
            return nothing
        end
        
        pkg_name = strip(parts[1])
        # Clean up the version string
        version = strip(replace(replace(parts[2], "\"" => ""), "~" => ""))
        
        return (pkg_name, version)
    catch e
        @warn "Failed to parse version line: $line" exception=e
        return nothing
    end
end

# Collect compatibility information with error tracking
compats = Dict{String,Tuple{Vector{String},Union{Nothing,Exception}}}()
for notebook in notebooks
    path = joinpath(@__DIR__, notebook)
    success, compat_lines, error = getcompat(path)
    if !success
        @warn "Failed to read compat info from $notebook" error
    end
    compats[notebook] = (compat_lines, error)
end

# Print detailed summary report
println("\n=== Summary Report ===")
println("Total notebooks found: ", length(notebooks))
println("Successfully updated: ", length(successful_updates))
println("Failed updates: ", length(failed_notebooks))

if !isempty(failed_notebooks)
    println("\nFailed Notebooks Details:")
    for (name, details) in failed_notebooks
        println("\n$name:")
        println("Error: ", details.error)
        println("Stack trace: ", first(details.trace))  # Show only first line of trace for brevity
    end
end

# Enhanced package version statistics
let
    println("\n=== Package Version Statistics ===")
    
    # Create a dictionary to store package versions by package name
    package_versions = Dict{String, Dict{String, Set{String}}}()
    
    # Collect all versions for each package
    for (notebook, (compat_lines, error)) in compats
        if !isnothing(error)
            continue
        end
        
        for line in compat_lines
            parsed = parse_version_line(line)
            isnothing(parsed) && continue
            
            pkg_name, version = parsed
            if !haskey(package_versions, pkg_name)
                package_versions[pkg_name] = Dict{String, Set{String}}()
            end
            
            # Get major version
            major = try
                split(version, ".")[1]
            catch
                "unknown"
            end
            
            if !haskey(package_versions[pkg_name], major)
                package_versions[pkg_name][major] = Set{String}()
            end
            push!(package_versions[pkg_name][major], notebook)
        end
    end
    
    # Print package statistics
    if isempty(package_versions)
        println("No package version information found in any notebook.")
    else
        for (pkg, versions) in package_versions
            println("\n$pkg:")
            for (major, notebooks) in versions
                println("  v$major.x.x used in $(length(notebooks)) notebooks:")
                for nb in notebooks
                    println("    - $nb")
                end
            end
        end
    end
end