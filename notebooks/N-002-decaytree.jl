### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 7cfb0c40-670b-11ef-3cee-81f5618ca090
begin
	using AbstractTrees
	using LinearAlgebra
	using OrderedCollections
end

# ╔═╡ a6ea3458-f7c8-467f-b5cc-3d7a0e7091f4
md"""
# Decay Tree in Julia
"""

# ╔═╡ 8e239522-ed97-4226-ae8f-59d9561715d5
t = ((1,2),(3,(4,5)));

# ╔═╡ c36eb901-b5dd-4ae6-9c23-294442f1299b
print_tree(t)

# ╔═╡ f6cfd8cd-5505-49bf-8d33-44440f63d0d0
map(PreOrderDFS(t)) do n
	ch = children(n)
	children_info = join(string.(collect(ch)), " and ")
	length(ch) != 0 && println("$n -> $children_info")
end

# ╔═╡ f2725edf-a6f5-4258-9517-2c6064a59408
PreOrderDFS(t) |> collect

# ╔═╡ 995b1880-391d-47c6-b9f3-e95068b5aaf5
md"""
## Example with four-vectors
"""

# ╔═╡ 3ff6b033-ffa9-42ac-8368-222ed4775932
momenta = [
	rand(4), # 1
	rand(4), # 1
	rand(4), # 1
	rand(4), # 1
	rand(4) # 1
]

# ╔═╡ 61c4e494-61e7-4ccd-8faa-c44219fff495
masssq(p) = p[4]^2-norm(p[1:3])^2

# ╔═╡ 76e806a3-6da7-4042-8a6f-3f45260ba4b1
masssq.(momenta)

# ╔═╡ 248e5da0-cb21-42c6-8e1f-55f2a5320131
function angles(tree, momenta)
	d = LittleDict()
	for n in PreOrderDFS(tree)
		all_children = children(n)
		length(all_children) == 0 && continue
		m1sq, m2sq = map(all_children) do child
			iv = collect(Leaves(child))
			p = sum(iv) do i
				momenta[i]
			end
			masssq(p)
		end
		m0sq = sum(i->momenta[i], collect(Leaves(n))) |> masssq
		d[n] = (; m1sq, m2sq, m0sq)
	end
	return d
end

# ╔═╡ 22359ecb-9307-4c77-88b6-9d12923fdd54
angles(t, momenta)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AbstractTrees = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OrderedCollections = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"

[compat]
AbstractTrees = "~0.4.5"
OrderedCollections = "~1.6.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "ed5c2af50587220f1f7c8b3d3bf900216999a252"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╟─a6ea3458-f7c8-467f-b5cc-3d7a0e7091f4
# ╠═7cfb0c40-670b-11ef-3cee-81f5618ca090
# ╠═8e239522-ed97-4226-ae8f-59d9561715d5
# ╠═c36eb901-b5dd-4ae6-9c23-294442f1299b
# ╠═f6cfd8cd-5505-49bf-8d33-44440f63d0d0
# ╠═f2725edf-a6f5-4258-9517-2c6064a59408
# ╟─995b1880-391d-47c6-b9f3-e95068b5aaf5
# ╠═3ff6b033-ffa9-42ac-8368-222ed4775932
# ╠═76e806a3-6da7-4042-8a6f-3f45260ba4b1
# ╠═61c4e494-61e7-4ccd-8faa-c44219fff495
# ╠═248e5da0-cb21-42c6-8e1f-55f2a5320131
# ╠═22359ecb-9307-4c77-88b6-9d12923fdd54
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
