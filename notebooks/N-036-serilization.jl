### A Pluto.jl notebook ###
# v0.19.29

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

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "b6f7799c28f0b99fac714bc997cf0c073ec48c06"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

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
