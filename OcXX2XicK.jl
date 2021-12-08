### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ b57fe6e0-96d9-11eb-1643-13cc8f32aa2b
begin
	using SymPy
	
	import PyCall
	PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	
	import_from(sympy.physics.wigner)
    import_from(sympy.physics.quantum.spin, (:WignerD,), typ=:Any)
end

# ╔═╡ b0a663a9-08cd-4334-8f9e-b4bb66e52145
using Plots

# ╔═╡ c98f5a1b-5f6b-4360-9a0a-ad5a0a5d1165
md"""
## Angular distribition in $\Omega_c^{**0}\to \Xi_c^{+}(\to pK^-\pi^+) K^-$
"""

# ╔═╡ f6a415a0-5630-49b9-a393-68bebd009b32
L(two_J,P) = div(two_J + (((mod(two_J+1,4)==2) == (P=='+')) ? 1 : -1), 2)

# ╔═╡ 5f5efe5c-c2c2-4e21-8556-fe229c2de192
two_JPs = [(1,'-'), (1,'+'), (3,'+'), (3,'-'), (5,'-')]

# ╔═╡ 4b7af024-d9af-49b4-9afa-12243c9f9491
[L(x...) for x in two_JPs]

# ╔═╡ de5507df-5b3d-4b74-8c5e-b1a58dfbef58
LJ(two_JP) = (L=L(two_JP...),J=Sym(two_JP[1])/2)

# ╔═╡ 419304a6-a6a1-4475-875b-23d55d0e42f0
md"""
### tests
"""

# ╔═╡ 2a9e6ec5-c65d-4bef-ab02-686a6aa0c910
md"""
### libraries and functions
"""

# ╔═╡ e1b15b5b-cfb3-435d-af99-a7efef32764f
theme(:wong, frame=:box, minorticks=true, grid=false)

# ╔═╡ 05bc24ae-a59a-452e-9ec1-3d5cfda7bf39
clg(j1,m1,j2,m2,j,m) = clebsch_gordan(Sym(j1),j2,j,m1,m2,m);

# ╔═╡ 4b2b2ea8-6b9e-4675-97e2-bcf92fcdf390
θ,z = @vars θ z positive=true

# ╔═╡ ecc592ce-0010-4188-bd8c-6fc83ed53128
dl(L,J,λ,ν) = sqrt(2L+Sym(1)) / sqrt(2J+Sym(1)) * sqrt(Sym(2)) *
	WignerD(Sym(1)/2,λ,ν,0,θ,0).doit() * clg(L,0,Sym(1)/2,λ,J,λ)

# ╔═╡ c13d9421-4997-45eb-959a-ef43ceaa16a7
I(L,J;α,β) = sympy.expand_trig(simplify(
		sum(abs2(dl(L,J,Sym(two_λ)/2,Sym(two_ν)/2))*(1+two_ν*α)/2*(1+two_λ*β)/2
			for two_λ in -1:2:1,
				two_ν in -1:2:1)))

# ╔═╡ c1597b63-b99a-408a-8339-c8d812dddcbb
Iz(two_JP;α=Sym(1),β=0.5) = I(LJ(two_JP)...;α=α,β=β).subs(sin(θ)^2, 1-z^2).subs(cos(θ)^2, z^2)

# ╔═╡ a2aeb283-4310-43ca-a319-121b87e89a29
Iz.(two_JPs)

# ╔═╡ Cell order:
# ╟─c98f5a1b-5f6b-4360-9a0a-ad5a0a5d1165
# ╠═ecc592ce-0010-4188-bd8c-6fc83ed53128
# ╠═c13d9421-4997-45eb-959a-ef43ceaa16a7
# ╠═f6a415a0-5630-49b9-a393-68bebd009b32
# ╠═5f5efe5c-c2c2-4e21-8556-fe229c2de192
# ╠═4b7af024-d9af-49b4-9afa-12243c9f9491
# ╠═de5507df-5b3d-4b74-8c5e-b1a58dfbef58
# ╠═c1597b63-b99a-408a-8339-c8d812dddcbb
# ╠═a2aeb283-4310-43ca-a319-121b87e89a29
# ╟─419304a6-a6a1-4475-875b-23d55d0e42f0
# ╟─2a9e6ec5-c65d-4bef-ab02-686a6aa0c910
# ╠═b57fe6e0-96d9-11eb-1643-13cc8f32aa2b
# ╠═b0a663a9-08cd-4334-8f9e-b4bb66e52145
# ╠═e1b15b5b-cfb3-435d-af99-a7efef32764f
# ╠═05bc24ae-a59a-452e-9ec1-3d5cfda7bf39
# ╠═4b2b2ea8-6b9e-4675-97e2-bcf92fcdf390
