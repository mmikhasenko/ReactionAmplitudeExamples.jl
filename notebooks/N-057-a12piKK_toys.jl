### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ d8489d70-cf9f-11ec-1f26-0b616e422d52
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"),
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("LaTeXStrings"),
		Pkg.PackageSpec("Parameters"),
		Pkg.PackageSpec("StaticArrays"),
		Pkg.PackageSpec("ThreadsX"),
		Pkg.PackageSpec("DelimitedFiles"),
	])
	# 
	using ThreeBodyDecay
	using Plots
	using LaTeXStrings
	using Parameters
	using StaticArrays
	using ThreadsX
	using DelimitedFiles
end

# ╔═╡ 47666396-b2dd-40c7-afe2-02fa25591963
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    lw=1, lab="", colorbar=false)

# ╔═╡ bb9a04d6-8fe9-4295-9a87-74200c485013
md"""
### Lineshape: standard Breit-Wigners
"""

# ╔═╡ cc71b82b-a1fa-43ff-879a-56aadda3b483
function BlattWeiss(z²,l)
	l==0 && return 1.0
	l==1 && return 1/(1+z²)
	l==2 && return 1/(9+3z²+z²^2)
	error("Not implemented, l=$l")
end

# ╔═╡ c135e8ca-3272-47bd-aad2-d7fe68b58352
md"""
### Specific case
"""

# ╔═╡ e357f507-0730-4a4b-9c5e-689adf4734be
begin
	const mπ = 0.14
	const mK = 0.497611
	const mKˣ, ΓKˣ = 0.892, 0.052
	const mf₀, Γf₀ = 0.996, 0.020
	const ma1 = 1.26
end ;

# ╔═╡ 11bd1e32-a909-4971-88c0-ca36d4f7942a
begin
	@with_kw struct RBW
		m::Float64
		mi::Float64
		mj::Float64
		Γ::Float64
		l::Int = 0
		R::Float64 = 1.5
	end
	(lineshape::RBW)(σ::Float64) = amplitude(lineshape, σ)
	# 
	import ThreeBodyDecay: amplitude
	function amplitude(lineshape::RBW, σ::Float64)
		@unpack m, Γ, mi, mj, l, R = lineshape
		sqrtσ = sqrt(σ)
		p = sqrtKallenFact(sqrtσ,mi,mj)/(2sqrtσ)
		p0 = sqrtKallenFact(m,mi,mj)/(2m)
		# 
		ff² = BlattWeiss((p*R)^2,l) / BlattWeiss((p0*R)^2,l)
		Γ_dep = Γ*ff²*(p/p0)^(2l+1)*m/sqrtσ
		return p^l*sqrt(ff²)*m*Γ/(m^2-σ-1im*m*Γ_dep)
	end
	@with_kw struct X3b{T}
		lineshape::T
		m0::Float64
		mk::Float64
		L::Int = 0
		r0::Float64 = 1.5
	end
	function (x::X3b)(σ::Float64)
		@unpack lineshape, m0, mk, r0, L = x
		sqrtσ = sqrt(σ)
		q = sqrtKallenFact(m0,sqrtσ,mk)/(2m0)
		ff² = BlattWeiss((q*r0)^2,L)
		return q^L*sqrt(ff²)*amplitude(lineshape, σ)
	end
	# 
	const Wave = NamedTuple{(:lineshape3b, :j0, :P0, :k, :lineshape, :j, :parity)}
	ms(m0) = ThreeBodyMasses(mπ,mK,mK; m0)
	# 
	function wave2decaychain(wave::Wave, s::Float64)
		@unpack lineshape3b, j0, P0 = wave
		@unpack k, lineshape, j, parity = wave
		# 
		m0 = sqrt(s)
		_ms = ms(m0)
		tbs = ThreeBodySystem(; ms=_ms, two_js=(0,0,0,j0) .|> x2)
		Ps = ('-','-','-',P0)
		d = DecayChainLS(2, lineshape; two_s=2j, parity, Ps, tbs)
		Xlineshape = X3b(; lineshape, m0, mk=_ms[k], L=div(d.HRk.two_ls[1],2))
		return DecayChainLS(k, Xlineshape; two_s=2j, parity, Ps, tbs)
	end
	function amplitude(wave::Wave, s::Float64, σs, λ::Int)
		@unpack lineshape3b = wave
		dX = wave2decaychain(wave, s)
		return Osym(dX,σs,λ)*amplitude(lineshape3b,s)
	end
end

# ╔═╡ 3e672231-c34c-4c99-8d68-2802607a20cc
function Osym(dc::T, σs, λ::Int) where T<:DecayChain
	dc.k == 1 && return amplitude(dc, σs, (0,0,0,2λ))
	dc.k == 2 && return amplitude(dc, σs, (0,0,0,2λ)) +
			minusone()^div(dc.two_s,2) * 
				amplitude(DecayChain(dc, k=3), σs, (0,0,0,2λ))
	error("Unaccounted case, k=$(dc.k)")
end

# ╔═╡ ceb3339d-6fbd-4643-8d46-da2abd6e1af7
begin
	BW_Kˣ = RBW(m=mKˣ, Γ=ΓKˣ, mi=mK, mj=mπ, l=1)
	BW_f₀ = RBW(m=mf₀, Γ=Γf₀, mi=mK, mj=mK, l=0)
end ;

# ╔═╡ 1aa2004e-0e5f-4a6e-9aaf-141d1ac51139
function intensity(x::X3b, σ::Float64)
	@unpack lineshape, m0, mk, mk = x
	@unpack mi, mj = lineshape
	sqrtσ = sqrt(σ)
	q = sqrtKallenFact(m0,sqrtσ,mk)/(2m0)
	p = sqrtKallenFact(sqrtσ,mi,mj)/(2sqrtσ)
	return abs2(lineshape(σ))*p*q/sqrtσ
end

# ╔═╡ 9b7ec391-7f41-40e7-b6ea-6176cc86fdfe
begin
	BW_a₁1260 = RBW(m=1.26, Γ=0.5, mi=2mπ, mj=mπ, l=0)
	BW_a₁1420 = RBW(m=1.42, Γ=0.15, mi=2mK, mj=mπ, l=1)
end ;

# ╔═╡ 67015bfb-b52f-4868-a4f4-c55904cfdde5
plot(ylab=L"|A|^2", xlab=L"m(KK\pi)\,\,(\mathrm{GeV})",
	plot(m->abs2(amplitude(BW_a₁1260,m^2)), 2mK, 2.5),
	plot(m->abs2(amplitude(BW_a₁1420,m^2)), 2mK, 2.5))

# ╔═╡ 01ea0651-7afc-40d1-95c5-5e4e91b9058e
begin
	struct model{N,T<:Wave}
		chains::StaticVector{N,T}
		couplings::StaticVector{N,ComplexF64}
	end
	model(dv::AbstractArray, cv::AbstractArray) =
		(N = length(dv); model(SVector{N}(dv), SVector{N}(cv .* (1+0im))))
end

# ╔═╡ 92c48924-1bb9-4eb4-a860-ccb80bf1dbcd
md"""
### An example
"""

# ╔═╡ ef92fad3-6242-45df-896a-306f30164e72
const waves = Wave[
	(lineshape3b=BW_a₁1260, j0=1, P0='+', k=2, lineshape=BW_Kˣ, j=1, parity='-')
	(lineshape3b=BW_a₁1420, j0=1, P0='+', k=1, lineshape=BW_f₀, j=0, parity='+')
] ;

# ╔═╡ e1fef9e1-bfcd-4b81-9247-ac467a3c3d3d
const defaultmodel = model(waves, [1.2,3.0*cis(π/2)])

# ╔═╡ 40473fac-21c4-4be7-8a3e-99fbede99b86
const jMax = maximum(getproperty.(waves, :j))

# ╔═╡ 88b7a1a1-f787-413f-a891-fa28f53b9378
intensity(m::model, s, σs) = sum(abs2, 
		sum(SVector([amplitude(d, s, σs, λ) for d in m.chains]) .* m.couplings)
	for λ in -jMax:jMax)

# ╔═╡ 93967929-62ef-44e5-8c89-d12f0505d8e3
plot(ylab=L"|A|^2", xlab=[L"m^2(KK)\,\,(\mathrm{GeV}^2)" L"m^2(K\pi)\,\,(\mathrm{GeV}^2)"],
	plot(σ->intensity(X3b(; lineshape=BW_f₀, m0=1.5, mk=mπ),σ), 4mK^2, 1.1),
	plot(σ->intensity(X3b(; lineshape=BW_Kˣ, m0=1.5, mk=mK),σ), (mπ+mK)^2, 1.1))

# ╔═╡ 89debcdf-8e45-4fd1-999c-0cfcde36a224
md"""
### Computation on a grid
"""

# ╔═╡ c4bc1803-12c7-443a-a14a-684a32e6c7a8
plot(ms(1.5), σs->intensity(defaultmodel, 1.5^2, σs), iσx=2, iσy=3)

# ╔═╡ a6f6f7f2-3ae2-4dc5-a742-01356c7d9390
md"""
### Numerical exercise
"""

# ╔═╡ 74bf268f-68fe-431c-9268-291c7ce792c6
begin
	torange(x,(x1,x2)) = x1+x*(x2-x1)
	function torange((x1,x2,x3); m0lims)
		m0 = torange(x1^2, m0lims)
		_ms = ms(m0)
		σ2 = torange(x2, lims2(_ms))
		σ3 = torange(x3, lims3(_ms))
		# 
		(; ms=_ms, σs=Invariants(_ms; σ2,σ3))
	end
end

# ╔═╡ ee4e1317-e391-4b7c-9b5c-41a470a76f96
data = filter!(
	mapslices(rand(3,100_000); dims=1) do x
		torange(x; m0lims = (2mK+mπ, 3.0))
	end[1,:]) do (ms, σs)
		inphrange(σs, ms)
end ;

# ╔═╡ ef5cfbb2-e086-4f1c-92ca-8e910e48dda4
length(data)

# ╔═╡ cfb646c1-3774-4fe4-8f74-afb861f75da7
Aiv = ThreadsX.collect(
    SVector([amplitude(d, ms.m0^2, σs, λ) for d in defaultmodel.chains])
    for λ in -jMax:jMax, (ms,σs) in data);

# ╔═╡ 6ca5bac5-fb4a-4c59-8686-4a19fc00159f
Iv = sum(Aiv, dims=1) do x
    abs2(sum(x .* defaultmodel.couplings))
end[1,:] ;

# ╔═╡ 28682328-e590-4070-ab15-fa866a8ca73f
begin
	Iv_i = []
	for i in 1:length(defaultmodel.couplings)
		cmap = [i==j for j in 1:length(defaultmodel.couplings)]
		_Iv = sum(Aiv, dims=1) do x
		    abs2(sum(x .* defaultmodel.couplings .* cmap))
		end[1,:] ;
		push!(Iv_i, _Iv) 
	end
end

# ╔═╡ 72b5debf-025b-4abd-9f8e-8fb9ea319515
fa11420 = sum(Iv_i[2]) / sum(Iv) ;

# ╔═╡ d4bd89ed-88ef-4967-b202-89847f68a087
md"""
#### Fraction of $a_1(1420)$ is $(round(fa11420*100; digits=1)) %
"""

# ╔═╡ 98dd240d-7bf7-4916-a511-e1bddfe79f22
let 
	bins = 100
	plot(xlab=L"m(KK\pi)", ylab=L"\# \mathrm{entries}")
	m0v = getproperty.(getindex.(data,:ms),:m0)
	stephist!(m0v; bins, weights=Iv, lab=L"\mathrm{total}")
	stephist!(m0v; bins, weights=Iv_i[1], lab=L"a_1(1260)")
	p = stephist!(m0v; bins, weights=Iv_i[2], lab=L"a_1(1420)")
	ymax = maximum(p[1][1][:y])
	lens!([1.1,1.8],[0, 2*ymax*fa11420], inset=(1,bbox(0.45,0.3,0.5,0.5)))
end

# ╔═╡ 971e1a90-aa4f-4c42-9c02-b6b3d3a58be8
md"""
## Save MC to disc
"""

# ╔═╡ 3f38c845-1fef-410a-bbb3-909e13428417
let
	m0v = getproperty.(getindex.(data, :ms), :m0)
	σ2v = getproperty.(getindex.(data, :σs), :σ2)
	σ3v = getproperty.(getindex.(data, :σs), :σ3)
	writedlm("a1toys.txt", [m0v σ2v σ3v Iv Iv_i[1] Iv_i[2]])
end

# ╔═╡ 7b5a5709-d482-459f-8bb8-a7a554eee4b0
md"""
#### File size: $(round(stat("a1toys.txt").size/2^20, digits=2)) Mb
"""

# ╔═╡ Cell order:
# ╠═d8489d70-cf9f-11ec-1f26-0b616e422d52
# ╠═47666396-b2dd-40c7-afe2-02fa25591963
# ╟─bb9a04d6-8fe9-4295-9a87-74200c485013
# ╠═cc71b82b-a1fa-43ff-879a-56aadda3b483
# ╠═11bd1e32-a909-4971-88c0-ca36d4f7942a
# ╠═3e672231-c34c-4c99-8d68-2802607a20cc
# ╟─c135e8ca-3272-47bd-aad2-d7fe68b58352
# ╠═e357f507-0730-4a4b-9c5e-689adf4734be
# ╠═ceb3339d-6fbd-4643-8d46-da2abd6e1af7
# ╠═1aa2004e-0e5f-4a6e-9aaf-141d1ac51139
# ╠═93967929-62ef-44e5-8c89-d12f0505d8e3
# ╠═9b7ec391-7f41-40e7-b6ea-6176cc86fdfe
# ╠═67015bfb-b52f-4868-a4f4-c55904cfdde5
# ╠═01ea0651-7afc-40d1-95c5-5e4e91b9058e
# ╟─92c48924-1bb9-4eb4-a860-ccb80bf1dbcd
# ╠═ef92fad3-6242-45df-896a-306f30164e72
# ╠═e1fef9e1-bfcd-4b81-9247-ac467a3c3d3d
# ╠═40473fac-21c4-4be7-8a3e-99fbede99b86
# ╠═88b7a1a1-f787-413f-a891-fa28f53b9378
# ╟─89debcdf-8e45-4fd1-999c-0cfcde36a224
# ╠═c4bc1803-12c7-443a-a14a-684a32e6c7a8
# ╟─a6f6f7f2-3ae2-4dc5-a742-01356c7d9390
# ╠═74bf268f-68fe-431c-9268-291c7ce792c6
# ╠═ee4e1317-e391-4b7c-9b5c-41a470a76f96
# ╠═ef5cfbb2-e086-4f1c-92ca-8e910e48dda4
# ╠═cfb646c1-3774-4fe4-8f74-afb861f75da7
# ╠═6ca5bac5-fb4a-4c59-8686-4a19fc00159f
# ╠═28682328-e590-4070-ab15-fa866a8ca73f
# ╠═72b5debf-025b-4abd-9f8e-8fb9ea319515
# ╟─d4bd89ed-88ef-4967-b202-89847f68a087
# ╠═98dd240d-7bf7-4916-a511-e1bddfe79f22
# ╟─971e1a90-aa4f-4c42-9c02-b6b3d3a58be8
# ╠═3f38c845-1fef-410a-bbb3-909e13428417
# ╟─7b5a5709-d482-459f-8bb8-a7a554eee4b0
