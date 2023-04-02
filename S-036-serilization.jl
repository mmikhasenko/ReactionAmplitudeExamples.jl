### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ e9502380-f2d2-11ec-0d1d-7bc10372ea64
using JSON

# ╔═╡ b27dad19-ffc0-4d53-a6fe-d6ba9b86d4d4
function writejson(path, obj)
    open(path, "w") do io
        JSON.print(io, obj)
    end
end

# ╔═╡ d38d861d-4720-49c3-a82c-a24eb2e5e0fd
N = 100

# ╔═╡ 3662d87a-2b0f-4332-a716-7f4a4e975831
xv = range(-1,1,length=N)

# ╔═╡ 40206149-e123-4810-896e-6b135858ef64
yv = range(-1,1,length=N)

# ╔═╡ 596434b3-e651-492f-b955-0cddcc450ef9
randzmatrix(N, fNaN=0.5) = map(x->x<fNaN ? NaN : x, rand(N,N))

# ╔═╡ f1293010-0d0a-478a-b696-6d8b10966441
filename = joinpath(@__DIR__,"gridded.json")

# ╔═╡ 0575d78e-f663-4e9f-be57-c9b77421a5c0
writejson(filename, Dict(
	:x => xv,
	:y => yv,
	:z1 => randzmatrix(N),
	:z2 => randzmatrix(N),
	:z3 => randzmatrix(N),
	:z0 => randzmatrix(N)
))

# ╔═╡ da260979-1533-45e2-9a93-6c7f9ef5336e
md"""
#### File size: $(round(stat(filename).size / 2^20, digits=2))Mb
"""

# ╔═╡ 507d6047-ae53-44df-b1ee-b0d57871f808
rm(filename)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"

[compat]
JSON = "~0.21.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═e9502380-f2d2-11ec-0d1d-7bc10372ea64
# ╠═b27dad19-ffc0-4d53-a6fe-d6ba9b86d4d4
# ╠═d38d861d-4720-49c3-a82c-a24eb2e5e0fd
# ╠═3662d87a-2b0f-4332-a716-7f4a4e975831
# ╠═40206149-e123-4810-896e-6b135858ef64
# ╠═596434b3-e651-492f-b955-0cddcc450ef9
# ╠═f1293010-0d0a-478a-b696-6d8b10966441
# ╠═0575d78e-f663-4e9f-be57-c9b77421a5c0
# ╟─da260979-1533-45e2-9a93-6c7f9ef5336e
# ╠═507d6047-ae53-44df-b1ee-b0d57871f808
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
