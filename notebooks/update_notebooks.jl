import Pluto
import Pkg

notebooks = readdir(@__DIR__)
filter!(notebooks) do name
    contains(name, "N-0")
end

map(notebooks) do name
    Pluto.activate_notebook_environment(joinpath("notebooks", name))
    Pkg.update()
end


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

let
    paths = joinpath.("notebooks", notebooks)
    compats = paths .=> getcompat.(paths)
end