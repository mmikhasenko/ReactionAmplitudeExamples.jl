#!/usr/bin/env julia
# validate_notebook_index.jl
#
# This script validates that all Pluto notebooks in the notebooks/ directory
# are properly linked in the README.md file.
#
# It checks for:
# 1. Notebooks linked in README.md but not found on disk
# 2. Notebooks on disk but not linked in README.md
# 3. Duplicate links in README.md
#
# README.md is considered the single source of truth.
# If any issues are found, it will exit with a non-zero status code.

function main()
    base_dir = dirname(dirname(abspath(@__FILE__)))
    readme_path = joinpath(base_dir, "README.md")
    notebooks_dir = joinpath(base_dir, "notebooks")

    # 1. Read README.md and extract notebook paths using regex
    readme_content = read(readme_path, String)
    linked_notebooks_relative = String[]
    # Regex to find markdown links like (notebooks/N-XXX-something.jl) or (notebooks/L-XXX-something.jl)
    regex = r"\(notebooks/([NL]-\d{3}-.*?\.jl)\)"
    for match in eachmatch(regex, readme_content)
        push!(linked_notebooks_relative, "notebooks/" * match.captures[1])
    end
    # Keep duplicates for now to check later

    # 2. Check if linked notebooks exist on disk
    nonexistent_files = filter(path -> !isfile(joinpath(base_dir, path)), linked_notebooks_relative) |> unique

    # 3. Get all actual .jl files in the notebooks directory
    actual_notebook_files = filter(f -> endswith(f, ".jl") && !endswith(f, ".plutostate"), readdir(notebooks_dir))
    actual_notebook_paths = ["notebooks/$file" for file in actual_notebook_files]

    # Skip certain files like test_notebooks.jl if necessary
    excluded_files = ["notebooks/test_notebooks.jl"]
    actual_notebook_paths = filter(path -> !(path in excluded_files), actual_notebook_paths)

    # 4. Check for notebooks on disk but not linked in README
    unlinked_notebooks = setdiff(actual_notebook_paths, linked_notebooks_relative)

    # 5. Check for duplicate links in README
    duplicate_links = filter(nb -> count(x -> x == nb, linked_notebooks_relative) > 1, linked_notebooks_relative) |> unique

    # Print summary and details
    issues_found = false

    if !isempty(nonexistent_files)
        issues_found = true
        println("\n❌ ERROR: The following notebooks are linked in README.md but do not exist on disk:")
        for nb in nonexistent_files
            println("  - $nb")
        end
        println("\nSuggestion: Remove these links from README.md or ensure the files exist.")
    end

    if !isempty(unlinked_notebooks)
        issues_found = true
        println("\n❌ ERROR: The following notebooks exist on disk but are not linked in README.md:")
        for nb in unlinked_notebooks
            println("  - $nb")
        end
        println("\nSuggestion: Add links for these notebooks to the appropriate section in README.md.")
    end

    if !isempty(duplicate_links)
        issues_found = true
        println("\n❌ ERROR: The following notebooks are linked more than once in README.md:")
        for nb in duplicate_links
            println("  - $nb")
        end
        println("\nSuggestion: Ensure each notebook is linked exactly once in README.md.")
    end

    if issues_found
        println("\n❌ Validation failed! Please fix the issues above.")
        exit(1)
    else
        println("\n✅ README.md notebook links are consistent with the notebooks directory.")
        println("   Found $(length(unique(linked_notebooks_relative))) unique notebook links in README.md.")
        println("   Found $(length(actual_notebook_paths)) notebooks on disk (excluding test files).")
        println("   All linked notebooks exist, and all existing notebooks are linked exactly once.")
        exit(0)
    end
end

# Run the script
main()