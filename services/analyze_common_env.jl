#!/usr/bin/env julia
"""
    analyze_common_env.jl

Analyzes all Pluto notebooks with PLUTO_PROJECT_TOML_CONTENTS and generates
a minimal common Project.toml that can run all notebooks.

Usage:
    julia services/analyze_common_env.jl [--output path/to/Project.toml]
"""

using Pkg
using TOML
using Glob

"""
    find_notebooks_with_project_toml(root_dir::String = ".")

Find all notebooks containing PLUTO_PROJECT_TOML_CONTENTS.
"""
function find_notebooks_with_project_toml(root_dir::String=".")
    notebooks = String[]

    # Find all .jl files in notebooks/ and literate/ directories
    notebook_dirs = [
        joinpath(root_dir, "notebooks"),
        joinpath(root_dir, "literate"),
        joinpath(root_dir, "jupyter")
    ]

    for dir in notebook_dirs
        if isdir(dir)
            # Find all .jl files recursively (including root level)
            for file in vcat(glob("*.jl", dir), glob("**/*.jl", dir))
                try
                    content = read(file, String)
                    # Check if it's a Pluto notebook and has PLUTO_PROJECT_TOML_CONTENTS
                    if startswith(content, "### A Pluto.jl notebook ###") &&
                       occursin("PLUTO_PROJECT_TOML_CONTENTS", content)
                        push!(notebooks, file)
                    end
                catch e
                    @warn "Could not read $file: $e"
                end
            end
        end
    end

    return notebooks
end

"""
    extract_project_toml_content(notebook_path::String)

Extract PLUTO_PROJECT_TOML_CONTENTS from a notebook.
"""
function extract_project_toml_content(notebook_path::String)
    try
        content = read(notebook_path, String)

        # Find PLUTO_PROJECT_TOML_CONTENTS = """..."""
        pattern = r"PLUTO_PROJECT_TOML_CONTENTS\s*=\s*\"\"\"(.*?)\"\"\""s
        match_result = match(pattern, content)

        if isnothing(match_result)
            return nothing
        end

        toml_content = String(strip(match_result.captures[1]))
        # Return nothing if content is empty
        return isempty(toml_content) ? nothing : toml_content
    catch e
        @warn "Error extracting Project.toml from $notebook_path: $e"
        return nothing
    end
end

"""
    parse_project_toml(toml_string::String)

Parse a TOML string into a dictionary.
"""
function parse_project_toml(toml_string::String)
    try
        return TOML.parse(toml_string)
    catch e
        @warn "Error parsing TOML: $e"
        return nothing
    end
end

"""
    merge_project_tomls(project_tomls::Vector{Dict})

Merge multiple Project.toml dictionaries into a common one.
"""
function merge_project_tomls(project_tomls::Vector{Dict})
    merged = Dict{String,Any}()

    # Collect all dependencies
    all_deps = Dict{String,String}()  # name => uuid

    # Collect all compat entries
    all_compat = Dict{String,Vector{String}}()  # name => [version specs]

    for proj in project_tomls
        # Merge dependencies
        if haskey(proj, "deps")
            for (name, uuid) in proj["deps"]
                # Use first UUID found (they should all be the same for same package)
                if !haskey(all_deps, name)
                    all_deps[name] = uuid
                elseif all_deps[name] != uuid
                    @warn "UUID mismatch for $name: $(all_deps[name]) vs $uuid"
                end
            end
        end

        # Merge compat entries
        if haskey(proj, "compat")
            for (name, version_spec) in proj["compat"]
                if !haskey(all_compat, name)
                    all_compat[name] = String[]
                end
                # Store version spec as string
                push!(all_compat[name], string(version_spec))
            end
        end
    end

    # Build merged Project.toml
    merged["deps"] = all_deps

    # For compat, we need to handle version ranges
    # This is simplified - in practice, you might want to use Pkg to resolve compatible versions
    if !isempty(all_compat)
        merged_compat = Dict{String,String}()
        for (name, version_specs) in all_compat
            # If all specs are the same, use that
            unique_specs = unique(version_specs)
            if length(unique_specs) == 1
                merged_compat[name] = unique_specs[1]
            else
                # Multiple different version specs - try to find common range
                # For now, we'll use the first one and warn
                # In practice, you'd want to use Pkg.Operations to resolve compatible versions
                @warn "Multiple version specs for $name: $unique_specs. Using first: $(unique_specs[1])"
                merged_compat[name] = unique_specs[1]
            end
        end
        merged["compat"] = merged_compat
    end

    return merged
end

"""
    generate_project_toml_string(merged::Dict)

Generate a Project.toml string from a dictionary.
"""
function generate_project_toml_string(merged::Dict)
    lines = String[]

    # Dependencies section
    push!(lines, "[deps]")
    deps = merged["deps"]
    for name in sort(collect(keys(deps)))
        uuid = deps[name]
        push!(lines, "$name = \"$uuid\"")
    end

    # Compat section (if present)
    if haskey(merged, "compat") && !isempty(merged["compat"])
        push!(lines, "")
        push!(lines, "[compat]")
        compat = merged["compat"]
        for name in sort(collect(keys(compat)))
            version_spec = compat[name]
            push!(lines, "$name = \"$version_spec\"")
        end
    end

    return join(lines, "\n") * "\n"
end

function main()
    # Parse command line arguments
    output_path = "common_Project.toml"
    if length(ARGS) > 0
        if ARGS[1] == "--output" && length(ARGS) > 1
            output_path = ARGS[2]
        elseif startswith(ARGS[1], "--")
            println("Usage: julia analyze_common_env.jl [--output path/to/Project.toml]")
            return
        else
            output_path = ARGS[1]
        end
    end

    # Change to project root directory
    project_root = dirname(dirname(abspath(@__FILE__)))
    cd(project_root)

    println("🔍 Finding notebooks with PLUTO_PROJECT_TOML_CONTENTS...")
    notebooks = find_notebooks_with_project_toml(".")
    println("   Found $(length(notebooks)) notebooks")

    println("\n📦 Extracting Project.toml contents...")
    project_tomls = Dict[]
    failed_notebooks = String[]

    for notebook in notebooks
        toml_content = extract_project_toml_content(notebook)
        if isnothing(toml_content)
            push!(failed_notebooks, notebook)
            continue
        end

        parsed = parse_project_toml(toml_content)
        if isnothing(parsed)
            push!(failed_notebooks, notebook)
            continue
        end

        push!(project_tomls, parsed)
        println("   ✓ $(basename(notebook))")
    end

    if !isempty(failed_notebooks)
        println("\n⚠️  Failed to extract Project.toml from:")
        for nb in failed_notebooks
            println("   - $nb")
        end
    end

    if isempty(project_tomls)
        error("No valid Project.toml contents found!")
    end

    println("\n🔀 Merging $(length(project_tomls)) Project.toml files...")
    merged = merge_project_tomls(project_tomls)

    println("\n📊 Summary before resolution:")
    println("   Total dependencies: $(length(merged["deps"]))")
    if haskey(merged, "compat")
        println("   Packages with compat entries: $(length(merged["compat"]))")
    end

    println("\n🔧 Resolving dependencies using Pkg...")
    # Use Pkg to resolve dependencies properly and generate Manifest
    temp_env = mktempdir()
    manifest_output_path = if endswith(output_path, "Project.toml")
        replace(output_path, "Project.toml" => "Manifest.toml")
    else
        replace(output_path, ".toml" => "_Manifest.toml")
    end
    manifest_content = nothing
    project_toml_string = generate_project_toml_string(merged)  # Initialize with merged version

    try
        Pkg.activate(temp_env)

        # Write merged Project.toml to temp environment
        temp_project = joinpath(temp_env, "Project.toml")
        project_toml_string = generate_project_toml_string(merged)
        write(temp_project, project_toml_string)

        # Let Pkg resolve all dependencies and generate Manifest
        Pkg.resolve()
        Pkg.instantiate()

        # Read back the resolved Project.toml
        resolved_project = TOML.parsefile(temp_project)

        # Read the generated Manifest.toml
        temp_manifest = joinpath(temp_env, "Manifest.toml")
        if isfile(temp_manifest)
            manifest_content = read(temp_manifest, String)
            println("   ✓ Dependencies resolved and Manifest generated")
        else
            @warn "Manifest.toml not generated"
        end

        # Generate final Project.toml string from resolved version
        project_toml_string = generate_project_toml_string(resolved_project)

        # Write Manifest if available
        if !isnothing(manifest_content)
            write(manifest_output_path, manifest_content)
            println("   ✓ Manifest written to: $manifest_output_path")
        end

    catch e
        @warn "Pkg resolution failed, using merged version: $e"
        # Fall back to merged version if resolution fails
        project_toml_string = generate_project_toml_string(merged)
    finally
        # Clean up temp environment
        rm(temp_env; recursive=true, force=true)
        # Restore original environment
        Pkg.activate(".")
    end

    println("\n💾 Generating common Project.toml...")

    # Write to file
    write(output_path, project_toml_string)
    println("   ✓ Written to: $output_path")

    # Also print to stdout
    println("\n" * "="^80)
    println("Generated Project.toml:")
    println("="^80)
    println(project_toml_string)

    println("\n✅ Done! Generated common environment files:")
    println("   - Project.toml: $output_path")
    if isfile(manifest_output_path)
        println("   - Manifest.toml: $manifest_output_path")
        println("\n   Both files should be used together when updating notebooks.")
    end
    println("   Note: Version compatibilities may need manual adjustment if conflicts exist.")
end

# Run if executed directly
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

