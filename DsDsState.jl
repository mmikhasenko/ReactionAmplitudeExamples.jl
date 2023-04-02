### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 3c395aaa-0660-4d61-bac2-aff30fe8768b
begin
	import Pkg
	Pkg.add(
	Pkg.PackageSpec(url="https://github.com/mmikhasenko/EffectiveRangeExpansion.jl"))
	Pkg.add("Measurements")
	Pkg.add("Plots")
	Pkg.add("Parameters")
	Pkg.add("LinearAlgebra")
	Pkg.add("Distributions")
	# 
	using Measurements
	using Parameters
	using Plots
	using LinearAlgebra
	using EffectiveRangeExpansion
	using Distributions
end ;

# ╔═╡ b2775c19-294a-45e4-a106-3bb41e32c8ea
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lw=1.2, lab="", colorbar=false)

# ╔═╡ c4f58fe6-b4c6-4e51-80d6-9dc7e2f145ca
begin
	const mDs⁺ = 1.96835 # GeV
	const mD⁺ = 1.86966 # GeV
end

# ╔═╡ 35021cc3-8f3c-4e5d-ad8d-975bc6b078c1
begin
	sqrtλ(a,b,c) = sqrt((a-(b+c)))*sqrt((a-(b-c)))*sqrt((a+(b+c)))*sqrt((a+(b-c)))
	q(m,m1,m2) = sqrtλ(m,m1,m2) / (2m)
	ρ(m,m1,m2) = 2*q(m,m1,m2) / m
	ρ1(m) = ρ(m,mDs⁺,mDs⁺)
	ρ2(m) = ρ(m,mD⁺,mD⁺)
end

# ╔═╡ 9b6a3cd5-1e0f-4c06-9725-bf0459e29ce9
begin
	e2m(e) = 2mDs⁺+e*1e-3
	m2e(m) = (m-2mDs⁺)*1e3
end

# ╔═╡ 59c76afe-fdd7-11ec-3cd4-53454db53f8d
function D(m; p)
	@unpack M0, g1, g2 = p
	return M0^2- m^2 - 1im*M0*(g1*ρ1(m)+g2*ρ2(m))
end

# ╔═╡ d312becb-1549-4e44-aecd-034ae417c024
function intensity(e; p)
	e < 0 && return 0.0
	m = e2m(e)
	ρ1(m)*m / abs2(D(m; p))
end

# ╔═╡ 217384a2-6697-4c86-857c-057bc26977c3
md"""
## Fit results and intensity
"""

# ╔═╡ 08143b1d-027d-481a-9173-665b620b86cf
begin
	M0_fit = 3.949 ± 0.017 # GeV
	g1_fit = 0.45 ± 1.22 # GeV
	g_fit2 = 0.14 ± 0.11 # GeV
end ;

# ╔═╡ 9f57ecaa-a050-4ff2-b09a-fe2cf4d4dd63
p0 = (M0_fit.val, g1_fit.val, g_fit2.val) |> NamedTuple{(:M0,:g1,:g2)}

# ╔═╡ f3456168-ebcf-4ebf-b343-c40e333e2071
begin
	plot(layout=grid(2,1, heights=(0.95,0.05)), link=:x)
	plot!(e->intensity(e; p=p0), 0, 200, xlab="Δm (MeV)", title="intensity spectrum")
	plot!(sp=2, yaxis=nothing, frame=:axes, xlab="m (GeV)",
		xticks = (0:50:200, round.(e2m.(0:50:200); digits=2)))
end

# ╔═╡ ed447ad4-6207-4184-bc1e-22cd6a774af5
md"""
## Effective range parameters
"""

# ╔═╡ 367607a5-2dcd-4851-8f3c-b2a3e29eb69e
begin
	targetf(Δe) = D(e2m(Δe); p=p0)
	k(Δe) = ρ1(e2m(Δe))*e2m(Δe)/2
end

# ╔═╡ bbc07cae-20f5-4aff-b065-8666db0424bc
const ere_method = ComplexBranchPointExpansion(CircularIntegral(0.1))

# ╔═╡ b0dd5a1a-a120-4fa6-b34b-76b08c6e845a
erp = let
	v = effectiverangeexpansion(targetf, k, ere_method)
	ERP(; NamedTuple{(:N, :a⁻¹, :r)}(v)...)
end

# ╔═╡ 6fc9980a-e878-41f4-a50c-a0df533e37dd
md"""
## Sampling parameters
"""

# ╔═╡ caa3c578-15b7-4424-a73a-9484b2bd5403
Corr = Symmetric(
	[   1.        0.95089658  -0.93626351
	0.95089658  1.          -0.97610506
	-0.93626351 -0.97610506  1.    ])

# ╔═╡ 933005c9-c2a6-45e1-a911-78213f4ab910
DG = begin
	μ = [M0_fit.val, g1_fit.val, g_fit2.val]
	σ = [M0_fit.err, g1_fit.err, g_fit2.err]
	Σ = σ' .* Corr .* σ
	MvNormal(μ, Σ)
end ;

# ╔═╡ 68d717fa-dbe2-40c4-a8cd-b9ddd7c685fc
pv = mapslices(rand(DG, 10_000); dims=1) do v
	NamedTuple{(:M0,:g1,:g2)}(v)
end[1,:] ;

# ╔═╡ a9c6881c-0890-465f-86f0-72dfac64f16a
let
	plot(title="Intensity", xlab="Δm (MeV)")
	for p in pv[1:100:end]
		xv = range(-20, 100, 300)
		yv = intensity.(xv; p)
		plot!(xv, yv ./ sum(yv), lc=1, lw=0.3)
	end
	plot!()
end

# ╔═╡ 0f9ef1b0-01b0-49ca-9c3e-75e333e88f01
erpv = [effectiverangeexpansion(Δe->D(e2m(Δe); p), k, ere_method) for p in pv] ;

# ╔═╡ 283f5839-9fd3-45e8-8927-0e54e3cd0cbf
begin
	scatter(title="scattering parameters", xlab="ℜ a⁻¹ (MeV)", ylab="ℜ r (MeV)",
		real.(getproperty.(erpv, :a⁻¹)),
		real.(getproperty.(erpv, :r)),
		xlim = (-0.3,0.3), ylim=(-10,10), ms=1)
	hline!([0.0], lab="", lc=1)
	vline!([0.0], lab="", lc=1)
end

# ╔═╡ eb63fd94-da91-43fa-a3df-4f20e1224ba7
BW3915(e; p = (m=3.919, Γ=13e-3)) = 1/(p.m^2-e^2-1im*p.m*p.Γ)

# ╔═╡ ce219e0c-7248-4c70-bb19-82645a2214fa
begin
	plot(layout=grid(2,1, heights=(0.95,0.05)), link=:x)
	plot!(e->abs2(BW3915(e2m(e)))*ρ1(e2m(e))*e2m(e), 0, 200, xlab="Δm (MeV)", title="intensity spectrum")
	plot!(sp=2, yaxis=nothing, frame=:axes, xlab="m (GeV)",
		xticks = (0:50:200, round.(e2m.(0:50:200); digits=2)))
end

# ╔═╡ Cell order:
# ╠═3c395aaa-0660-4d61-bac2-aff30fe8768b
# ╠═b2775c19-294a-45e4-a106-3bb41e32c8ea
# ╠═c4f58fe6-b4c6-4e51-80d6-9dc7e2f145ca
# ╠═35021cc3-8f3c-4e5d-ad8d-975bc6b078c1
# ╠═9b6a3cd5-1e0f-4c06-9725-bf0459e29ce9
# ╠═59c76afe-fdd7-11ec-3cd4-53454db53f8d
# ╠═d312becb-1549-4e44-aecd-034ae417c024
# ╟─217384a2-6697-4c86-857c-057bc26977c3
# ╠═08143b1d-027d-481a-9173-665b620b86cf
# ╠═9f57ecaa-a050-4ff2-b09a-fe2cf4d4dd63
# ╠═f3456168-ebcf-4ebf-b343-c40e333e2071
# ╟─ed447ad4-6207-4184-bc1e-22cd6a774af5
# ╠═367607a5-2dcd-4851-8f3c-b2a3e29eb69e
# ╠═bbc07cae-20f5-4aff-b065-8666db0424bc
# ╠═b0dd5a1a-a120-4fa6-b34b-76b08c6e845a
# ╟─6fc9980a-e878-41f4-a50c-a0df533e37dd
# ╠═caa3c578-15b7-4424-a73a-9484b2bd5403
# ╠═933005c9-c2a6-45e1-a911-78213f4ab910
# ╠═68d717fa-dbe2-40c4-a8cd-b9ddd7c685fc
# ╠═a9c6881c-0890-465f-86f0-72dfac64f16a
# ╠═0f9ef1b0-01b0-49ca-9c3e-75e333e88f01
# ╠═283f5839-9fd3-45e8-8927-0e54e3cd0cbf
# ╠═eb63fd94-da91-43fa-a3df-4f20e1224ba7
# ╠═ce219e0c-7248-4c70-bb19-82645a2214fa
