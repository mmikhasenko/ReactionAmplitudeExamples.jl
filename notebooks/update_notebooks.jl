import Pluto
import Pkg
import StatsBase: countmap


notebooks = readdir(@__DIR__)
filter!(notebooks) do name
    contains(name, r"N-0")
end

map(notebooks) do name
    println("\n\n⋆⋆⋆ ", name, " ⋆⋆⋆\n\n")
    Pluto.activate_notebook_environment(joinpath("notebooks", name))
    if Pkg.API.Context().env.manifest.manifest_format != v"2.0.0"
        println("🐺🐺🐺 Upgrading ", name)
        Pkg.upgrade_manifest()
        Pkg.resolve()
    end
    Pkg.update()
end
Pkg.activate(joinpath(@__DIR__, ".."))

function getcompat(file)
    compat = []
    start_fill = false
    for line in readlines(file)
        if line == "[compat]"
            start_fill = true
        end
        if start_fill && line == "\"\"\""
            start_fill = false
        end
        start_fill && push!(compat, line)
    end
    compat[2:end]
end

compats = let
    paths = joinpath.("notebooks", notebooks)
    paths .=> getcompat.(paths)
end

# We look over all required packages and
# make sure that the major vertion is the same between the notebooks.
let
    allversion = getindex.(compats, 2)
    allversion_combined = sort(vcat(allversion...))
    cut_to_major = map(allversion_combined) do p
        parts = split(p, " = ")
        parts[1] * " = " * parts[2][1:3]
    end
    remove_repetitions = Set(cut_to_major) |> collect
    packages_distinct_major = getindex.(split.(remove_repetitions, " = "), 1)
    repetitions = countmap(packages_distinct_major)
    for p in repetitions
        p[2] > 1 && @warn "Distinct major versions for $(p[1]) 😨"
    end
end