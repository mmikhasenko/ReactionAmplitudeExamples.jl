### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ f1634740-cfe1-11ed-3804-b74bfa23da7a
using SymPy

# ╔═╡ afc4e979-e5f0-495b-a10e-c1bca0845dc7
md"""
# Cartesian indices for $3\pi$ decays

We work out covariant amplitudes for $X^i \to \pi^j \pi^k \pi^k$ decays for arbitrary isospin using the cartesian basis. The amplitude is related to the charge basis by linear transformation.

We focus on the $J^{PC} = 2^{++}$ sector.
"""

# ╔═╡ 21ec2e15-a984-4816-bbb1-08efbb1016ee
const sqrt2 = sqrt(Sym(2)) ;

# ╔═╡ 8d2f9797-2593-45c6-ace0-cda0b1afe808
begin
	M = sympy.IndexedBase("\\mathcal{M}") # total matrix  element
	@syms i j k l # cartesian indices
	# 
	δ(i,j) = sympy.KroneckerDelta(i,j)
	# 
	Mv = @syms M_s M_t M_u
	# 
	p = sympy.IndexedBase("p")  # 4-vector
	Bv = @syms B_s B_t B_u  # single variable amplitudes
end;

# ╔═╡ 70fdc5c8-a8b2-492d-b9d7-f98043f0bd9a
amplitude_cartesian = Mv[1]*δ(i,j)*δ(k,l) + Mv[2]*δ(i,k)*δ(j,l) + Mv[3]*δ(i,l)*δ(k,j)

# ╔═╡ 7596ff8c-2778-4a05-ad1e-ac40cc52cc3d
md"""
Define pion struct for convenience. Essentually, we only need the charge.
"""

# ╔═╡ 7a1ec453-4743-40bb-b9ec-f4982cd66ef1
begin
	struct Pion
		charge::Int
	end
	Pion(f::(typeof(-))) = Pion(-1)
	Pion(f::(typeof(+))) = Pion(1)
	up(charge) = charge==0 ? '⁰' : (charge==1 ? '⁺' : '⁻')
	Base.show(io::IO, pi::Pion) = print(io, Symbol("π"*up(pi.charge)))
	Base.conj(pi::Pion) = Pion(-pi.charge)
end

# ╔═╡ 660a389c-7e77-4480-a0ab-42254292dc57
md"""
The reaction is a decay $X\to 3\pi$. For the mapping of the charges to the cartesian coordinates, the final-state fields are conjugated. Instead of conjugating coefficients, I conjugate the pion charges.
"""

# ╔═╡ 9a4835f2-a86e-4970-bb0c-0b5ed41a870a
begin
	const Reaction = Pair{Pion, Tuple{Pion, Pion, Pion}}
	charges(r::Reaction) = (r[1], conj.(r[2])...)	
end ;

# ╔═╡ c0bfa54c-8206-4ba7-a7a9-4b00fbdb2d11
function project(pi,c)
	# 
	X = [1/sqrt2 0   1/sqrt2
	     -1im/sqrt2 0  1im/sqrt2
	     0       0         1] .|> Sym
	# 
	X[pi.charge+2, c]
end

# ╔═╡ adda2d2b-69c3-4d5a-b795-06d0df493ddf
function charge2cartesian(expr, r::Reaction)
	ch = charges(r)
	sum(Iterators.product(Iterators.repeated(1:3,4)...)) do indices in 
		coef = prod(project.(ch, indices))
		coef * expr.xreplace(Dict((i,j,k,l) .=> indices))
	end
end

# ╔═╡ 8239f355-7b1b-40c9-85d9-f85b121408f6
const reaction_tu = Pion(-) => (Pion(+),Pion(-),Pion(-))

# ╔═╡ 16950a25-dfde-4984-b794-bb818c04b5da
charge2cartesian(M[i,j,k,l], reaction_tu)

# ╔═╡ ce919c17-d3f5-45f1-8647-8060ff462b9b
md"""
The model for the matrix element is

$M_s(s,t,u) = (p_2+p_3) B(s,t,u) + (p_2-p_3) C(s,t,u)\,$

where the scalar functions $B$ and $C$ are expressed via a single-varible function B(s)

```math
\begin{align*}
B(s,t,u) = B(t)-B(u)\,,\\
C(s,t,u) = B(t)+B(u)\,.
\end{align*}
```

The amplitudes $M_t(s,t,u)$, and $M_u(s,t,u)$ are defined by the circular permutations $(s\to t \to u \to s)$ applied to the function $M_s(s,t,u)$.
"""

# ╔═╡ 0df233e5-135d-4cf1-84b4-d4b210438564
transform2BC = map(1:3) do k
	i,j = mod.(k.+[0,1],3) .+ 1
	# 
	Mv[k] => (p[Sym(i)]+p[Sym(j)]) * (Bv[i]-Bv[j]) +
		(p[Sym(i)]-p[Sym(j)]) * (Bv[i]+Bv[j]) |> simplify
end |> Dict

# ╔═╡ a9c7d255-d13e-40c7-8a77-4214c98d01d8
map((Pion(-) => (Pion(-),Pion(-),Pion(+)),
	Pion(-) => (Pion(+),Pion(-),Pion(-)),
	Pion(-) => (Pion(-),Pion(+),Pion(-)))) do r
	
	(r, charge2cartesian(amplitude_cartesian, r).xreplace(transform2BC))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
SymPy = "~1.1.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0-rc1"
manifest_format = "2.0"
project_hash = "fc0c2fe4421f3fe529dfc9943e1bbac9ca0a667b"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CommonEq]]
git-tree-sha1 = "d1beba82ceee6dc0fce8cb6b80bf600bbde66381"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.0"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "e32a90da027ca45d84678b826fffd3110bb3fc90"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.8.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "62f417f6ad727987c755549e9cd88c46578da562"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.95.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.SymPy]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "PyCall", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "fcb24df16e451cfa8e1e1217edfd92054f75d49d"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "1.1.8"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.4.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─afc4e979-e5f0-495b-a10e-c1bca0845dc7
# ╠═f1634740-cfe1-11ed-3804-b74bfa23da7a
# ╠═21ec2e15-a984-4816-bbb1-08efbb1016ee
# ╠═8d2f9797-2593-45c6-ace0-cda0b1afe808
# ╠═70fdc5c8-a8b2-492d-b9d7-f98043f0bd9a
# ╟─7596ff8c-2778-4a05-ad1e-ac40cc52cc3d
# ╠═7a1ec453-4743-40bb-b9ec-f4982cd66ef1
# ╟─660a389c-7e77-4480-a0ab-42254292dc57
# ╠═9a4835f2-a86e-4970-bb0c-0b5ed41a870a
# ╠═c0bfa54c-8206-4ba7-a7a9-4b00fbdb2d11
# ╠═adda2d2b-69c3-4d5a-b795-06d0df493ddf
# ╠═8239f355-7b1b-40c9-85d9-f85b121408f6
# ╠═16950a25-dfde-4984-b794-bb818c04b5da
# ╟─ce919c17-d3f5-45f1-8647-8060ff462b9b
# ╠═0df233e5-135d-4cf1-84b4-d4b210438564
# ╠═a9c7d255-d13e-40c7-8a77-4214c98d01d8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
