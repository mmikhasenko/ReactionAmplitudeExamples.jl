
files = readdir()

fs = filter(files) do f
    f[end-2:end] == ".jl" && f != "summary.jl" && isfile(f)
end

ns = filter(fs) do f
    first(eachline(f)) == "### A Pluto.jl notebook ###"
end

for n in ns
    newname = "N" * n[2:end]
    # println(n * " -> " * newname)
    mv(n, newname)
end

