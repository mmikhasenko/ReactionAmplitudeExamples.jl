### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ d8489d70-cf9f-11ec-1f26-0b616e422d52
begin
	import Pkg
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"),
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("LaTeXStrings"),
		Pkg.PackageSpec("Parameters"),
		Pkg.PackageSpec("StaticArrays"),
		Pkg.PackageSpec("ThreadsX"),
	])
	# 
	using ThreeBodyDecay
	using Plots
	using LaTeXStrings
	using Parameters
	using StaticArrays
	using ThreadsX
end

# ╔═╡ 47666396-b2dd-40c7-afe2-02fa25591963
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    lw=1, lab="", colorbar=false,
	xlab=L"m^2\,\,(\mathrm{GeV}^2)",
	ylab=L"m^2\,\,(\mathrm{GeV}^2)") # _{\pi^+ \pi^-}

# ╔═╡ cc71b82b-a1fa-43ff-879a-56aadda3b483
function BlattWeiss(z²,l)
	l==0 && return 1.0
	l==1 && return 1/(1+z²)
	l==2 && return 1/(9+3z²+z²^2)
	error("Not implemented, l=$l")
end

# ╔═╡ bb9a04d6-8fe9-4295-9a87-74200c485013
md"""
### Lineshape: standard Breit-Wigners
"""

# ╔═╡ 11bd1e32-a909-4971-88c0-ca36d4f7942a
begin
	@with_kw struct RBW
		m::Float64
		mi::Float64
		mj::Float64
		Γ::Float64
		l::Int = 0.0
		R::Float64 = 1.5
	end
	(lineshape::RBW)(σ::Float64) = amplitude(lineshape, σ)
	# 
	import ThreeBodyDecay: amplitude
	function amplitude(lineshape::RBW, σ::Float64)
		@unpack m, mi, mj, Γ, l, R = lineshape
		sqrtσ = sqrt(σ)
		p = sqrtKallenFact(sqrtσ,mi,mj)/(2sqrtσ)
		p0 = sqrtKallenFact(m,mi,mj)/(2m)
		# 
		R = 1.5 # 1/GeV
		ff² = BlattWeiss((p*R)^2,l) / BlattWeiss((p0*R)^2,l)
		Γ = Γ*ff²*(p/p0)^(2l+1)*m/sqrtσ
		p*sqrt(ff²)/(m^2-σ-1im*m*Γ)
	end
end

# ╔═╡ c1333000-1441-41e7-8085-e388a434c389
sqrtKallenFact(5,1,1), sqrt(Kallen(5^2,1,1))

# ╔═╡ c135e8ca-3272-47bd-aad2-d7fe68b58352
md"""
### Specific case
"""

# ╔═╡ ef3117f7-8198-4016-88f0-49ba6a8dd240
0.497611*2

# ╔═╡ e357f507-0730-4a4b-9c5e-689adf4734be
begin
	const mπ = 0.14
	const mK = 0.497611
	const mKˣ, ΓKˣ = 0.892, 0.052
	const mf₀, Γf₀ = 0.996, 0.010
	const ma1 = 1.46
end ;

# ╔═╡ 21922eab-95bc-436e-8b1d-f532cdfa99e8
const ms = ThreeBodyMasses(mπ,mK,mK; m0=ma1)

# ╔═╡ 9a9c76c4-2055-4f71-bb40-d96bc0cb82d8
const j0 = 1 #

# ╔═╡ fc6deed2-3dd4-4b82-9ea7-52d3ffd3af1e
const tbs_a1 = ThreeBodySystem(; ms, two_js=(0,0,0,j0) .|> x2)

# ╔═╡ 89831e2a-245f-47d8-9f4a-464e6c292cce
const Ps_a1 = ('-','-','-','+')

# ╔═╡ ceb3339d-6fbd-4643-8d46-da2abd6e1af7
begin
	BW_Kˣ = RBW(m=mKˣ, Γ=ΓKˣ, mi=mK, mj=mπ, l=1)
	BW_f₀ = RBW(m=mf₀, Γ=Γf₀, mi=mK, mj=mK, l=0)
end ;

# ╔═╡ 93967929-62ef-44e5-8c89-d12f0505d8e3
plot(ylab="|A|²",
	plot(σ->abs2(amplitude(BW_f₀,σ)), lims1(ms)...),
	plot(σ->abs2(amplitude(BW_Kˣ,σ)), lims2(ms)...))

# ╔═╡ f3bd334d-d6e6-4a52-b76c-fba25511fe96
const dcv = [
	DecayChainLS(2, BW_Kˣ; two_s=1|>x2, parity='-', Ps=Ps_a1, tbs=tbs_a1),
	DecayChainLS(1, BW_f₀; two_s=0|>x2, parity='+', Ps=Ps_a1, tbs=tbs_a1)]

# ╔═╡ 01ea0651-7afc-40d1-95c5-5e4e91b9058e
begin
	struct model{N,T<:DecayChain}
		chains::StaticVector{N,T}
		couplings::StaticVector{N,ComplexF64}
	end
	model(dv::AbstractArray, cv::AbstractArray) =
		(N = length(dv); model(SVector{N}(dv), SVector{N}(cv .* (1+0im))))
end

# ╔═╡ e1fef9e1-bfcd-4b81-9247-ac467a3c3d3d
const defaultmodel = model(dcv, [2.0,-10.0im])

# ╔═╡ 3e672231-c34c-4c99-8d68-2802607a20cc
function Osym(dc::T, σs, λ::Int) where T<:DecayChain
	dc.k == 1 && return amplitude(dc, σs, (0,0,0,2λ))
	dc.k == 2 && return amplitude(dc, σs, (0,0,0,2λ)) -
			minusone()^λ * amplitude(DecayChain(dc, k=3), σs, (0,0,0,2λ))
	error("Unaccounted case, k=$(dc.k)")
end

# ╔═╡ 88b7a1a1-f787-413f-a891-fa28f53b9378
intensity(m::model, σs) = sum(abs2, 
		sum(SVector([Osym(d, σs, λ) for d in m.chains]) .* m.couplings)
	for λ in -j0:j0)

# ╔═╡ 89debcdf-8e45-4fd1-999c-0cfcde36a224
md"""
### Computation on a grid
"""

# ╔═╡ c4bc1803-12c7-443a-a14a-684a32e6c7a8
plot(ms, Base.Fix1(intensity, defaultmodel), iσx=2, iσy=3)

# ╔═╡ a6f6f7f2-3ae2-4dc5-a742-01356c7d9390
md"""
### Numerical exercise
"""

# ╔═╡ 30fb7bea-ea06-4f9d-8d54-b6060ce049cd
const data = flatDalitzPlotSample(ms; Nev=500_000) ;

# ╔═╡ cfb646c1-3774-4fe4-8f74-afb861f75da7
Aiv = ThreadsX.collect(
    SVector([Osym(d, σs, λ) for d in defaultmodel.chains])
    for λ in -j0:j0, σs in data);

# ╔═╡ 6ca5bac5-fb4a-4c59-8686-4a19fc00159f
Iv = sum(Aiv, dims=1) do x
    abs2(sum(x .* defaultmodel.couplings))
end[1,:] ;

# ╔═╡ 9438ffc3-537a-4efc-b407-d5e26393855d
I0 = ThreadsX.sum(Aiv) do x
    abs2(sum(x .* defaultmodel.couplings))
end

# ╔═╡ 3cc87fac-5e88-47bf-b15f-0758275b9ab7
histogram2d(
	getproperty.(data,:σ2),
	getproperty.(data,:σ3),
	weights=Iv, bins=40)

# ╔═╡ Cell order:
# ╠═d8489d70-cf9f-11ec-1f26-0b616e422d52
# ╠═47666396-b2dd-40c7-afe2-02fa25591963
# ╠═cc71b82b-a1fa-43ff-879a-56aadda3b483
# ╟─bb9a04d6-8fe9-4295-9a87-74200c485013
# ╠═11bd1e32-a909-4971-88c0-ca36d4f7942a
# ╠═c1333000-1441-41e7-8085-e388a434c389
# ╟─c135e8ca-3272-47bd-aad2-d7fe68b58352
# ╠═ef3117f7-8198-4016-88f0-49ba6a8dd240
# ╠═e357f507-0730-4a4b-9c5e-689adf4734be
# ╠═21922eab-95bc-436e-8b1d-f532cdfa99e8
# ╠═9a9c76c4-2055-4f71-bb40-d96bc0cb82d8
# ╠═fc6deed2-3dd4-4b82-9ea7-52d3ffd3af1e
# ╠═89831e2a-245f-47d8-9f4a-464e6c292cce
# ╠═ceb3339d-6fbd-4643-8d46-da2abd6e1af7
# ╠═93967929-62ef-44e5-8c89-d12f0505d8e3
# ╠═f3bd334d-d6e6-4a52-b76c-fba25511fe96
# ╠═01ea0651-7afc-40d1-95c5-5e4e91b9058e
# ╠═e1fef9e1-bfcd-4b81-9247-ac467a3c3d3d
# ╠═3e672231-c34c-4c99-8d68-2802607a20cc
# ╠═88b7a1a1-f787-413f-a891-fa28f53b9378
# ╟─89debcdf-8e45-4fd1-999c-0cfcde36a224
# ╠═c4bc1803-12c7-443a-a14a-684a32e6c7a8
# ╟─a6f6f7f2-3ae2-4dc5-a742-01356c7d9390
# ╠═30fb7bea-ea06-4f9d-8d54-b6060ce049cd
# ╠═cfb646c1-3774-4fe4-8f74-afb861f75da7
# ╠═6ca5bac5-fb4a-4c59-8686-4a19fc00159f
# ╠═9438ffc3-537a-4efc-b407-d5e26393855d
# ╠═3cc87fac-5e88-47bf-b15f-0758275b9ab7
