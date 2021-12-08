### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ a230fcb2-98a4-11eb-0ad0-c964bf424978
begin
	using SymPy
	
	import PyCall
	PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	
	import_from(sympy.physics.wigner)
    import_from(sympy.physics.quantum.spin, (:WignerD,), typ=:Any)
end

# ╔═╡ 5e2fd427-5f64-4c1e-b96a-5d5db9e433fa
begin
	using LinearAlgebra
	using Plots
end

# ╔═╡ cbcbdbcd-1f3d-4d9d-ac05-dfc4533e1d9e
md"""
# Amplitude for $e^+ e^- \to J/\psi \to \pi^+\pi^-\pi^0$
"""

# ╔═╡ 1103298c-e2fa-48b4-8bd9-daa5211f3e61
md"""
$A = \sum_{\nu} D_{1,\nu}^{1*}(0,\beta,\phi) O^{\nu}(\sigma_1,\sigma_3)$
"""

# ╔═╡ 614a6c9b-bc36-40e6-bead-8bfd66ff3a8e
β, ϕ, z = @vars β ϕ z positive=true

# ╔═╡ 3fbc1152-5d12-4d97-9578-192f1de5af66
a, b, c = @vars a b c

# ╔═╡ b966b6d9-b4dd-4e23-aadc-27b4e592bfb3
wignerd(j,m1,m2,θ) = WignerD(j,m1,m2,0,θ,0).doit()

# ╔═╡ ca073733-333d-49b3-8635-7aceb1b48a68
A = sympy.trigsimp(dot([WignerD(1,1,ν,0,β,ϕ).doit() for ν in -1:1], [1, 0, 1]))

# ╔═╡ f8d4e5fe-8a78-42cb-bc7c-831f846032b2
expand(A)

# ╔═╡ 5dde48cf-c7c0-4a09-a382-1f42908ce070
f = lambdify(A.subs(cos(β),z).subs(sin(β),sqrt(1-z^2)))

# ╔═╡ a6d6e3ad-9f37-4d0c-98f1-e1bb470aa3e9
let
	xv = range(-1,1, length=100)
	yv = range(-π,π, length=110)
	calv = abs2.(f.(xv',yv))
	heatmap(xv, yv, calv, c=:viridis, 
		xlab="cosβ", ylab="ϕ")
end

# ╔═╡ e2a257ce-7a48-4051-9464-6899e4e72e86
cg(j1,m1,j2,m2,j,m) = clebsch_gordan(Sym(j1),j2,j,m1,m2,m)

# ╔═╡ 7e08d73a-5bb0-4945-8ff6-888e128f870b
md"""
```math
\int_{-\pi}^{\pi} (\cos^2(\phi) + \cos^2(\theta) \sin^2(\phi) )\frac{\mathrm{d}\cos\theta}{2} = 
\cos^2(\phi) + \frac{1}{3}\sin^2(\phi)
```
"""

# ╔═╡ 257ab71f-b43c-464b-870c-41e7acb72470
Iϕ(ϕ) = cos(ϕ)^2+sin(ϕ)^2 / 3

# ╔═╡ e175a5f6-c89b-42af-842b-7510d1541511
plot(Iϕ, -π, π, ylim=(0,Inf))

# ╔═╡ be318ca8-a16e-4238-81fe-3c7af9398cf4
md"""
## Dalitz Plot distribution
$O^\nu = \sum_{k=1} d_{\nu,\lambda}^1(\theta_{k(1)})
	H_{\lambda}\,
	d_{\lambda,0}(\theta_{ij})$
where we use the circular notations, the chain-1 is the reference chain: $\hat{\theta}_{1(1)} = 0$
"""

# ╔═╡ c570ef7e-388e-4edc-81c6-4ee50dabdc90
θ, θhat = @vars θ θhat positive=true;

# ╔═╡ 7a854543-3ff6-461e-b2f2-6ccbb4f6a08f
cross_term(j,ν=1) = sympy.expand_trig(simplify(
	sum(wignerd(1,ν,λ,θhat)*cg(j,0,j,λ,1,λ)*wignerd(j,λ,0,θ) for λ in -1:1)))

# ╔═╡ 8188dd5c-c62a-4d3c-a864-59360dc32a43
[cross_term(1,ν) for ν in -1:1]

# ╔═╡ 65158ef9-86a6-4e94-8566-e02282cf2e28
simplify.([cross_term(3,ν) for ν in -1:1])

# ╔═╡ 6783d495-d1c0-4f72-ab44-74f2dbf4fb0f
md"""
## Simplification of the polynomials
"""

# ╔═╡ bd01f05a-ea8d-4978-af4b-7b26235152f6
s, σ1, σ2, σ3 = @vars s σ_1 σ_2 σ_3 positive=true;

# ╔═╡ b0433a76-93a5-4f2a-a644-1865d5c07042
s, σ1, σ2, σ3

# ╔═╡ 38c5cdc0-e084-4805-ab60-f7e76f5f50de
Φ = σ1*σ2*σ3 - (s-1)^2

# ╔═╡ 4911bcb9-ed74-4ae8-8cdf-b2d28d7cbdc4
λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x

# ╔═╡ 9f5d02e3-c8e8-413d-812c-afb0f64bcece
simplify(-5*σ1*Φ.subs(σ2,3+s-σ1-σ3) + λ(σ1,1,1)*λ(s,σ1,1))/σ1

# ╔═╡ Cell order:
# ╟─cbcbdbcd-1f3d-4d9d-ac05-dfc4533e1d9e
# ╠═a230fcb2-98a4-11eb-0ad0-c964bf424978
# ╠═1103298c-e2fa-48b4-8bd9-daa5211f3e61
# ╠═614a6c9b-bc36-40e6-bead-8bfd66ff3a8e
# ╠═3fbc1152-5d12-4d97-9578-192f1de5af66
# ╠═5e2fd427-5f64-4c1e-b96a-5d5db9e433fa
# ╠═b966b6d9-b4dd-4e23-aadc-27b4e592bfb3
# ╠═ca073733-333d-49b3-8635-7aceb1b48a68
# ╠═f8d4e5fe-8a78-42cb-bc7c-831f846032b2
# ╠═5dde48cf-c7c0-4a09-a382-1f42908ce070
# ╠═a6d6e3ad-9f37-4d0c-98f1-e1bb470aa3e9
# ╠═e2a257ce-7a48-4051-9464-6899e4e72e86
# ╟─7e08d73a-5bb0-4945-8ff6-888e128f870b
# ╠═257ab71f-b43c-464b-870c-41e7acb72470
# ╠═e175a5f6-c89b-42af-842b-7510d1541511
# ╟─be318ca8-a16e-4238-81fe-3c7af9398cf4
# ╠═c570ef7e-388e-4edc-81c6-4ee50dabdc90
# ╠═7a854543-3ff6-461e-b2f2-6ccbb4f6a08f
# ╠═8188dd5c-c62a-4d3c-a864-59360dc32a43
# ╠═65158ef9-86a6-4e94-8566-e02282cf2e28
# ╟─6783d495-d1c0-4f72-ab44-74f2dbf4fb0f
# ╠═bd01f05a-ea8d-4978-af4b-7b26235152f6
# ╠═b0433a76-93a5-4f2a-a644-1865d5c07042
# ╠═38c5cdc0-e084-4805-ab60-f7e76f5f50de
# ╠═4911bcb9-ed74-4ae8-8cdf-b2d28d7cbdc4
# ╠═9f5d02e3-c8e8-413d-812c-afb0f64bcece
