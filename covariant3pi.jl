### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 70065006-585e-4d3d-982f-11eb4746f51a
begin
	using Symbolics
	using Combinatorics
	using LinearAlgebra
	using WignerSymbols
	# 
	using Symbolics.Latexify
end

# ╔═╡ 9268d453-b7f2-43f6-9af3-0ff54b658598
md"""
# Projection of covariant amplitudes

[arXiv 2212.11767](https://arxiv.org/abs/2212.11767)

In this notebook, we will calculate projection kernels related to the quantification of rescattering effects in three-pion systems. We focus on two decay processes with total $J^{PC} = 1^{-+}$, and $2^{++}$, which decay predominantly as $\rho\pi$ states.

By solving the Khuri--Treiman integral equations, we aim to investigate the significance of rescattering effects beyond two-body resonances and determine the minimum number of events required to unambiguously find these effects in future Dalitz-plot analyses.

This analysis will help us identify kinematic effects that either enhance or dilute the rescattering for the selected set of quantum numbers and various masses.

The notebook is in Julia. I am experimenting with the `Symbolics` package.
"""

# ╔═╡ 19e2509a-f7a1-49b4-b79f-aae345e29c8e
@variables g[1:4,1:4] KroneckerDelta[1:4, 1:4] ϵ[1:4, 1:4, 1:4, 1:4]

# ╔═╡ b6fbaf6d-1bda-4bc0-b963-f6072de4e401
begin
	evaluateϵ = Dict(
		Symbolics.scalarize(ϵ .=> [levicivita(vcat(indices...))
		for indices in Iterators.product(Symbolics.shape(ϵ)...)]))
	# 
	evaluateδ = Dict(Symbolics.scalarize(KroneckerDelta .=> [(i==j) for (i,j) in
		Iterators.product(Symbolics.shape(KroneckerDelta)...)]))
	# 
	evaluateg = Dict(
	Symbolics.scalarize(g .=> Diagonal([-1,-1,-1,1]))) ;
end ;

# ╔═╡ d9ddee24-58dc-4db5-a65d-6b17c40cd8ef
md"""
## Define variables
"""

# ╔═╡ c3a5cc00-c4a9-493e-818a-6be87ab5742c
md"""
tensors
"""

# ╔═╡ a8479a8f-5268-43be-8d0f-6f29a83fbc59
begin
	#tensors
	@variables ε[1:4, 1:4] p_1[1:4]  p_2[1:4] p_3[1:4] p_0[1:4] 
	# constants: masses, Mandelstam variable $s$, scattering angle $\theta$
	@variables s m_1 m_2 m_3 m_0 θ
	@variables k p  # three-momenta $k$, and $p$
	@variables j  # imaginary part $j$
end ;

# ╔═╡ 725245d1-ac4a-40c4-a17d-0f01f3137fa6
md"""
## Evaluate four-momenta
"""

# ╔═╡ 00af5e32-761f-4bcf-82f8-9b584958955d
evaluatep2 = Symbolics.scalarize(
	p_2 .=> [k*sin(θ), 0, k*cos(θ), (s+m_2^2-m_3^2) / 2 / sqrt(s)]) ;

# ╔═╡ 7e2aa8a0-a19a-42fb-8c59-e7520515cd02
evaluatep3 = Symbolics.scalarize(
	p_3 .=> [-k*sin(θ), 0, -k*cos(θ), (s-m_2^2+m_3^2) / 2 / sqrt(s)]) ;

# ╔═╡ 61f44078-8b6a-4307-bfb5-bd498c7d028b
evaluatep1 = Symbolics.scalarize(
	p_1 .=> [0, 0, -p, (m_0^2-s-m_1^2) / 2 / sqrt(s)]) ;

# ╔═╡ a27ada9f-a392-49be-a0a7-daa10b51c4b4
E = (m_0^2+s-m_1^2) / 2 / sqrt(s) ;

# ╔═╡ a592d167-58ed-4cab-85bd-73bf10b6a7f4
evaluatep0 = Symbolics.scalarize(
	p_0 .=> [0, 0, -p, E]) ;

# ╔═╡ 555f9b43-9de9-4ece-95f8-47ca4c978e81
tosymbols = Dict(evaluatep1..., evaluatep2..., evaluatep3..., evaluatep0...) ;

# ╔═╡ 1d23967b-95d6-4470-b04f-236f025879c4
md"""
Check that $p_1+p_2+p_3=p_0$
"""

# ╔═╡ a16cd38f-2be0-4864-8704-793fa9a668fa
p_1+p_2+p_3-p_0 |> Symbolics.scalarize .|> x->substitute(x, tosymbols) |> simplify

# ╔═╡ f8a7b92f-caca-4feb-9816-b344a90ef17c
md"""
## Spin-1 vector

The polarization covariant vector depends on helicity, $\lambda \in \{-1,0,1\}$

$\epsilon^\nu(\lambda)$

"""

# ╔═╡ 497861b9-738a-466c-868c-21b5c6f16d89
@variables ε_1[1:4,1:3];

# ╔═╡ 661db8e3-6014-42be-b320-36abbc6f626a
evaluateε1 = let
	hm1 = [+1,j,0,0] ./ Symbolics.Term(sqrt,2)
	hz = [0,0,E,-p] ./ m_0
	hp1 = [-1,j,0,0] ./ Symbolics.Term(sqrt,2)
	
	Dict(Symbolics.scalarize(ε_1 .=> [hm1 hz hp1]))
end ;

# ╔═╡ 95da1b0f-0633-4f81-8425-375a6d7bb7b8
extrarules = let
	r1 = @rule sqrt(~x)^2 => ~x
	r2 = @rule ((~x)/(~y))^2 => (~x)^2 / (~y)^2
	RuleSet([r1,r2])
end ;

# ╔═╡ 1a76919a-9d19-488e-889f-e9a5c07630a1
md"""
## Compute the amplitude
"""

# ╔═╡ 87b38591-46d8-45a1-80c2-1abf4d9a6685
md"""
$K_{\nu} = \epsilon_{\alpha\beta\mu\nu} p_1^\beta p_2^\mu  p_3^\nu$
"""

# ╔═╡ bb237c06-6918-4e25-8f4c-154ced42e02b
K = (Symbolics.@arrayop (α,) ϵ[α,β,μ,ν]*p_1[β]*p_2[μ]*p_3[ν] term=:K);

# ╔═╡ 9943d7bd-4d50-4e51-a6e4-cb913db6500f
A = (Symbolics.@arrayop (λ,) ε_1[μ,λ]*K[μ]) ;

# ╔═╡ 3db65143-6d6a-4b60-ac49-eb76a16dd59a
function formulate_equations(d::Union{Dict, Vector{Pair{A,B}} where {A,B}})
	t = """
	```math
	\\begin{align}
	"""
	for (k,v) in d
		t *= latexify(k; env=:raw)
		t *= "&= " 
		t *= latexify(v; env=:raw)
		t *= " \\\\"
	end
	t *= """
	\\end{align}
	```
	"""
	Markdown.parse(t)
end

# ╔═╡ 2a8152f9-da8b-466f-8a58-3cea831ec759
formulate_equations(evaluatep2)

# ╔═╡ 22927bd3-e88d-4b4d-b5c3-ad4b54b1acb2
formulate_equations(evaluatep3)

# ╔═╡ a28fe6d3-90c9-4698-ba0b-87b905ec1c3b
formulate_equations(evaluatep1)

# ╔═╡ 381f83a7-25f4-4bc9-8320-1dce81a24b8d
formulate_equations(evaluatep0)

# ╔═╡ 1685198d-4961-4561-aa87-3b8316c89075
["$(ε_1.value.name)(" .* string.([-1 0 1]) .*")" =>
		(ε_1 |> Symbolics.scalarize .|>
	x->substitute(x, evaluateε1))
] |> formulate_equations

# ╔═╡ 037426e9-a38a-4956-86c5-90d1b4a004c8
[("$(K.term)[$(K.output_idx[1])]") => K.expr] |> formulate_equations

# ╔═╡ 089e2368-606e-40d2-96bd-6cd4a6705f9f
[("A[$(A.output_idx[1])]") => A.expr] |> formulate_equations

# ╔═╡ d943f094-ed03-440f-a1fa-ba606a5bba62
A_sc = A |> Symbolics.scalarize .|>
	x->substitute(x, evaluateϵ) .|>
	x->substitute(x, evaluateε1) .|>
	x->substitute(x, tosymbols; fold=false) .|>
	simplify;

# ╔═╡ de32d0d3-5617-4c43-a96f-5cbd0298abb6
("A[" .* string.(-1:1) .*"]" .=> A_sc) |> formulate_equations

# ╔═╡ 1bfcd6f1-74e2-436c-b467-e836ed22da0d
md"""
## Spin-2 tensor

The polarization covariant vector depends on helicity, $\lambda \in \{-2,-1,0,1,2\}$.

It is computed using the spin-1 vectors and Clebsch-Gordan coefficients

$\epsilon^{\nu\tau}(\lambda) = \sum_{\lambda_1,\lambda_2} C_{\lambda_1,\lambda_2}^{1,1,2}\epsilon^{\tau}(\lambda_1)\epsilon^{\tau}(\lambda_2)$

"""

# ╔═╡ cb30d7bb-4989-4874-87c2-5981dbe3954a
@variables CG[1:3,1:3,1:5];

# ╔═╡ 86e53a4e-1956-41d5-81b7-e7b40cdebe1a
function convert2term(x::WignerSymbols.RationalRoots.RationalRoot)
	x == 0 && return Num(0)
	Symbolics.Term(sqrt, abs(x.signedsquare)) * sign(x)
end ;

# ╔═╡ 2eb39b39-9902-487a-9726-98066d42bafa
evaluateCG = Dict(
		Symbolics.scalarize(CG) .=> 
			[(clebschgordan(1,λ1-2,1,λ2-2,2,λ-3) |> convert2term) for (λ1,λ2,λ) in
		Iterators.product(Symbolics.shape(CG)...)]
) ;

# ╔═╡ 4d21e3ff-6fcb-4f21-82e4-5ee162bbf363
ε_2 = Symbolics.@arrayop (μ,ν,λ) ε_1[μ,λ_1] * ε_1[ν,λ_2] * CG[λ_1,λ_2,λ] term=:ε_2;

# ╔═╡ 0344bcfe-f0ff-45da-b339-2cc0b30a0339
[("$(ε_2.term)[$(ε_2.output_idx...),]") => ε_2.expr] |> formulate_equations

# ╔═╡ 1ea0ab96-2361-40d3-97b2-73c7984b4dad
ε2_sc = ε_2[:,:,1] |> Symbolics.scalarize .|>
	x->substitute(x, evaluateCG) .|> 
	x->substitute(x, evaluateε1) .|> 
	x->substitute(x, j^2 => -1) ;

# ╔═╡ 68304160-1524-42cd-ac97-f65158f65c2d
md"""
## The $J^{PC}=2^{++}$ sector

The matrix element reads

```math
A = \varepsilon^{\nu\tau}(\lambda)\, K_\nu\,
\big[
	B(s,t,u)(p_{2}+p_{3})_\tau + C(s,t,u)(p_{2}-p_{3})_\tau
\big]
```

The matrix element in the paper is defined in extra numerical factor, $\mathcal{M} = i\sqrt{2}A$

The expression is evaluated term by term as follows:

```math
A = B(s,t,u) A_2 + C(s,t,u) A_{2\prime}
```
where $A_2$ and $A_{2\prime}$ are given below.

```math
\begin{align*}
A_2(\lambda) &=  \varepsilon^{\nu\tau}(\lambda)\, K_\nu\, (p_{2}+p_{3})_\tau\,,\\
A_{2\prime}(\lambda) &=  \varepsilon^{\nu\tau}(\lambda)\, K_\nu\, (p_{2}-p_{3})_\tau\,.
\end{align*}
```
"""

# ╔═╡ c52957e2-45fe-492d-b086-a1f5b422344e
A_2 = Symbolics.@arrayop (λ,) ε_2[μ,ν,λ]*K[μ]*g[ν,ν′]*(p_2[ν′]+p_3[ν′]);

# ╔═╡ d087fcd9-15d3-4e33-9480-f88c2b0aa8c7
["A_2[$(A_2.output_idx[1])]" => A_2.expr] |> formulate_equations

# ╔═╡ 5d4cebc4-a156-461b-aa61-4334abbdec03
md"""
Shame on Symbolics.jl, I have to do manual simplification
"""

# ╔═╡ bbacb260-37a7-4314-be0b-dbaa64905cb7
sqrtunit = Symbolics.Term(sqrt,1//2) * Symbolics.Term(sqrt,2) ;

# ╔═╡ ba40ce15-4153-4976-a4b8-0106436888f3
A_2_sc = A_2 |> Symbolics.scalarize .|> 
	x->substitute(x, j^2 => -1) .|> 
	x->substitute(x, evaluateϵ) .|> 
	x->substitute(x, evaluateg) .|> 
	x->substitute(x, evaluateε1) .|> 
	x->substitute(x, evaluateCG; fold=false) .|> 
	x->substitute(x, tosymbols; fold=false) .|> 
	simplify .|>
	x->simplify(x / sqrtunit, extrarules);

# ╔═╡ 7d91316b-3197-4fba-bd4e-3104d5f191cf
("A_2[".*string.(-2:2).*"]" .=> A_2_sc) |> formulate_equations

# ╔═╡ 0622242c-1646-424a-b6fd-e981cab77ad0
md"""
### Cross channel
"""

# ╔═╡ 6e9f5dbd-fb4e-4930-9d9b-5eedbc70b997
A_2′ = Symbolics.@arrayop (λ,) ε_2[μ,ν,λ]*K[μ]*g[ν,ν′]*(p_1[ν′]+p_2[ν′]);

# ╔═╡ b871d212-8c1a-4caa-9411-fa35624f47b2
["A_2′[$(A_2′.output_idx[1])]" => A_2′.expr] |> formulate_equations

# ╔═╡ 9b1a2f13-0751-4c5c-af1f-8d0acc64c525
A_2′_sc = A_2′ |> Symbolics.scalarize .|> 
	x->substitute(x, j^2 => -1) .|> 
	x->substitute(x, evaluateϵ) .|> 
	x->substitute(x, evaluateg) .|> 
	x->substitute(x, evaluateε1) .|> 
	x->substitute(x, evaluateCG; fold=false) .|> 
	x->substitute(x, tosymbols; fold=false) .|> 
	simplify .|>
	x->simplify(x, extrarules) ;

# ╔═╡ 77fcee8a-7038-419e-95e4-c792af03f86f
# ("A_2′[".*string.(-2:2).*"]" .=> A_2′_sc) |> formulate_equations

# ╔═╡ 5b2b8f46-d1ad-4570-a05a-7cb44d9a8dff
md"""
Shape the equations nicely:
 
a) helicity-0 amplitude
"""

# ╔═╡ 4e63dcb4-d53d-4bcb-8708-58148bd3e795
["A_2′[0]" => A_2′_sc[3]] |> formulate_equations

# ╔═╡ e87daf8c-4f3c-4844-bfca-6ca77b804ac6
md"""
b) helicity-1 amplitude
"""

# ╔═╡ 3f572ed0-60d6-4242-8bf4-3f85f53e6bca
function specificshape1(amplitude)
	expr = simplify(amplitude / sqrtunit, extrarules)
	# 
	factor = j*p*k / (4m_0)
	unfac = (expr / factor ) |> simplify
	
	expr1 = k*sin(2θ)/2
	expr2 = p*sin(θ)
	# 
	xslit = sum(unfac.val.dict) do (k,v)
		hassin2θ(x::Symbolics.Mul) = sin(2θ) in keys(x.dict)
		k * v .* (hassin2θ(k) ? [1/expr1, 0] : [0,1/expr2])
	end .* [expr1, expr2] |> sum
	# 
	return xslit * factor
	
end

# ╔═╡ 9a78dbba-5e73-4c74-9ec0-395c32d734d9
begin
	xs = specificshape1.(A_2′_sc[[2,4]])
	("A_2′[".*string.([-1,1]).*"]" .=> xs) |> formulate_equations
end

# ╔═╡ b06aeb16-7a65-46c0-88d2-351dc0393c19
md"""
Type in latex:
```math
\frac{i p^2 k s}{4m_0} \left[ - \frac{s+m_3^2-m_2^2}{s}\sin\theta + \frac{s+m_0^2-m_1^2}{s} \frac{k}{p}\sin\theta \cos\theta \right]
```
"""

# ╔═╡ 442b7ecb-7795-4f09-b961-20213ce66a65
md"""
c) helicity-2 amplitude
"""

# ╔═╡ f54fb79c-faf0-4ff6-9765-6320e407b6dc
let
	xs = A_2′_sc[[1,5]]
	xs .*= sqrt(s)*sqrt(s)/s / Symbolics.Term(sqrt, 1//1)
	("A_2′[".*string.([-2,2]).*"]" .=> xs) |> formulate_equations
end

# ╔═╡ 7797b445-a17e-4e2c-8415-6bd42a0e69cb
md"""
Type in latex:
```math
-\frac{i p k^2 \sqrt{s}}{2} \sin^2\theta
```
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
WignerSymbols = "9f57e263-0b3d-5e2e-b1be-24f2bb48858b"

[compat]
Combinatorics = "~1.0.2"
Symbolics = "~4.13.0"
WignerSymbols = "~2.0.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0-rc1"
manifest_format = "2.0"
project_hash = "3975ded72b4a77797790ec147b6f1d88f55d7825"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "7772df04fda9bc25a44c9ef61e9dc7c92bb35d86"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.27.7"

[[deps.AbstractTrees]]
git-tree-sha1 = "52b3b436f8f73133d7bc3a6c71ee7ed6ab2ab754"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.3"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "6d0918cb9c0d3db7fe56bea2bc8638fc4014ac35"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.24"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "badccc4459ffffb6bce5628461119b7057dec32c"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.27"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "f12dc65aef03d0a49650b20b2fdaf184928fd886"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.5"

[[deps.ArrayInterfaceStaticArraysCore]]
deps = ["Adapt", "ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "93c8ba53d8d26e124a5a8d4ec914c3a16e6a0970"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.3"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "7fe6d92c4f281cf4ca6f2fba0ce7b299742da7ca"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.37"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.DataAPI]]
git-tree-sha1 = "46d2680e618f8abd007bce0c3026cb0c4a8f2032"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.12.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "a7756d098cbabec6b3ac44f369f74915e8cfd70a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.79"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "988e2db482abeb69efc76ae8b6eba2e93805ee70"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.15"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "d0fa82f39c2a5cdb3ee385ad52bc05c42cb4b9f0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.5"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "9a0472ec2f5409db243160a8b030f94c380167a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.6"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "a69dd6db8a809f78846ff259298678f0d6212180"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.34"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "a5e6e7f12607e90d71b09e6ce2c965e41b337968"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.1"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "6872f5ec8fd1a38880f027a26739d42dcda6691f"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "47f0f03eddecd7ad59c42b1dd46d5f42916aff63"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.11"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HalfIntegers]]
git-tree-sha1 = "00db638039558e6396b93e2702862d6a884ac50e"
uuid = "f0d1745a-41c9-11e9-1dd9-e5d34d218721"
version = "1.4.3"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "3f91cd3f56ea48d4d2a75c2a65455c5fc74fa347"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.LRUCache]]
git-tree-sha1 = "d862633ef6097461037a00a13f709a62ae4bdfdd"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.4.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "dae002226b59701dbafd7e2dd757df1bd83442fd"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.12.5"

[[deps.LambertW]]
git-tree-sha1 = "2d9f4009c486ef676646bca06419ac02061c088e"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.5"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

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
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

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

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0f39bc7f71abdff12ead4fc4a7d998fb2f3c171f"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.5"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "4d5917a26ca33c66c8e5ca3247bd163624d35493"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "393fc4d82a73c6fe0e2963dd7c882b09257be537"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "aa532179d4a643d4bd9f328589ca01fa20a0d197"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.1.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

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
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ForwardDiff"]
git-tree-sha1 = "2c88339bcfa011089f7538f89c5be780d0b558bb"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.6"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RationalRoots]]
git-tree-sha1 = "52315cf3098691c1416a356925027af5ab5bf548"
uuid = "308eb6b3-cc68-5ff3-9e97-c3c4da4fa681"
version = "0.2.0"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "18c35ed630d7229c5584b945641a73ca83fb5213"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.2"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "ZygoteRules"]
git-tree-sha1 = "53d040e68b5afff59a9eb1f692d9a08369b07f61"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.34.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "50314d2ef65fce648975a8e80ae6d8409ebbf835"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.5"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "StaticArraysCore", "Statistics", "Tables"]
git-tree-sha1 = "0f016d69ed6df4ec438e468986036a75493c3261"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.79.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "c35b107b61e7f34fa3f124026f2a9be97dea9e1c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.3"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "ffc098086f35909741f71ce21d03dadf0d2bfa76"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.11"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "6b764c160547240d868be4e961a5037f47ad7379"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "027b43d312f6d52187bb16c2d4f0588ddb8c4bb2"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.11"

[[deps.Symbolics]]
deps = ["ArrayInterfaceCore", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "718328e81b547ef86ebf56fbc8716e6caea17b00"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.13.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "2d7164f7b8a066bcfa6224e67736ce0eb54aef5b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.9.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "34e6bcf36b9ed5d56489600cf9f3c16843fa2aa2"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.11"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f2fd3f288dfc6f507b0c3a2eb3bac009251e548b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.22"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c42fa452a60f022e9e087823b47e5a5f8adc53d5"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.75"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WignerSymbols]]
deps = ["HalfIntegers", "LRUCache", "Primes", "RationalRoots"]
git-tree-sha1 = "960e5f708871c1d9a28a7f1dbcaf4e0ee34ee960"
uuid = "9f57e263-0b3d-5e2e-b1be-24f2bb48858b"
version = "2.0.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

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
# ╟─9268d453-b7f2-43f6-9af3-0ff54b658598
# ╠═70065006-585e-4d3d-982f-11eb4746f51a
# ╠═19e2509a-f7a1-49b4-b79f-aae345e29c8e
# ╠═b6fbaf6d-1bda-4bc0-b963-f6072de4e401
# ╠═3db65143-6d6a-4b60-ac49-eb76a16dd59a
# ╟─d9ddee24-58dc-4db5-a65d-6b17c40cd8ef
# ╟─c3a5cc00-c4a9-493e-818a-6be87ab5742c
# ╠═a8479a8f-5268-43be-8d0f-6f29a83fbc59
# ╟─725245d1-ac4a-40c4-a17d-0f01f3137fa6
# ╠═00af5e32-761f-4bcf-82f8-9b584958955d
# ╠═2a8152f9-da8b-466f-8a58-3cea831ec759
# ╠═7e2aa8a0-a19a-42fb-8c59-e7520515cd02
# ╠═22927bd3-e88d-4b4d-b5c3-ad4b54b1acb2
# ╠═61f44078-8b6a-4307-bfb5-bd498c7d028b
# ╠═a28fe6d3-90c9-4698-ba0b-87b905ec1c3b
# ╠═a27ada9f-a392-49be-a0a7-daa10b51c4b4
# ╠═a592d167-58ed-4cab-85bd-73bf10b6a7f4
# ╠═381f83a7-25f4-4bc9-8320-1dce81a24b8d
# ╠═555f9b43-9de9-4ece-95f8-47ca4c978e81
# ╟─1d23967b-95d6-4470-b04f-236f025879c4
# ╠═a16cd38f-2be0-4864-8704-793fa9a668fa
# ╟─f8a7b92f-caca-4feb-9816-b344a90ef17c
# ╠═497861b9-738a-466c-868c-21b5c6f16d89
# ╠═661db8e3-6014-42be-b320-36abbc6f626a
# ╠═1685198d-4961-4561-aa87-3b8316c89075
# ╠═95da1b0f-0633-4f81-8425-375a6d7bb7b8
# ╟─1a76919a-9d19-488e-889f-e9a5c07630a1
# ╠═87b38591-46d8-45a1-80c2-1abf4d9a6685
# ╠═bb237c06-6918-4e25-8f4c-154ced42e02b
# ╠═037426e9-a38a-4956-86c5-90d1b4a004c8
# ╠═9943d7bd-4d50-4e51-a6e4-cb913db6500f
# ╠═089e2368-606e-40d2-96bd-6cd4a6705f9f
# ╠═d943f094-ed03-440f-a1fa-ba606a5bba62
# ╠═de32d0d3-5617-4c43-a96f-5cbd0298abb6
# ╟─1bfcd6f1-74e2-436c-b467-e836ed22da0d
# ╠═cb30d7bb-4989-4874-87c2-5981dbe3954a
# ╠═86e53a4e-1956-41d5-81b7-e7b40cdebe1a
# ╠═2eb39b39-9902-487a-9726-98066d42bafa
# ╠═4d21e3ff-6fcb-4f21-82e4-5ee162bbf363
# ╠═0344bcfe-f0ff-45da-b339-2cc0b30a0339
# ╠═1ea0ab96-2361-40d3-97b2-73c7984b4dad
# ╟─68304160-1524-42cd-ac97-f65158f65c2d
# ╠═c52957e2-45fe-492d-b086-a1f5b422344e
# ╠═d087fcd9-15d3-4e33-9480-f88c2b0aa8c7
# ╟─5d4cebc4-a156-461b-aa61-4334abbdec03
# ╠═bbacb260-37a7-4314-be0b-dbaa64905cb7
# ╠═ba40ce15-4153-4976-a4b8-0106436888f3
# ╠═7d91316b-3197-4fba-bd4e-3104d5f191cf
# ╟─0622242c-1646-424a-b6fd-e981cab77ad0
# ╠═6e9f5dbd-fb4e-4930-9d9b-5eedbc70b997
# ╟─b871d212-8c1a-4caa-9411-fa35624f47b2
# ╠═9b1a2f13-0751-4c5c-af1f-8d0acc64c525
# ╠═77fcee8a-7038-419e-95e4-c792af03f86f
# ╟─5b2b8f46-d1ad-4570-a05a-7cb44d9a8dff
# ╠═4e63dcb4-d53d-4bcb-8708-58148bd3e795
# ╟─e87daf8c-4f3c-4844-bfca-6ca77b804ac6
# ╠═3f572ed0-60d6-4242-8bf4-3f85f53e6bca
# ╠═9a78dbba-5e73-4c74-9ec0-395c32d734d9
# ╟─b06aeb16-7a65-46c0-88d2-351dc0393c19
# ╟─442b7ecb-7795-4f09-b961-20213ce66a65
# ╠═f54fb79c-faf0-4ff6-9765-6320e407b6dc
# ╟─7797b445-a17e-4e2c-8415-6bd42a0e69cb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
