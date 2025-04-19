### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a4f34cb4-170b-11f0-1df4-114a305c6853
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/NumericalDistributions.jl"),
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("ComponentArrays"),
		Pkg.PackageSpec("FHist")
		])
	using NumericalDistributions
	using NumericalDistributions.Distributions
	# 
	using Plots
	using FHist
	using ComponentArrays
end

# ╔═╡ 647abad9-5891-4154-8f91-2163b5825434
theme(:boxed, fillalpha=0.3)

# ╔═╡ 29e9dda3-6b55-4dc7-af7e-f2d0b20bd3ae
g = NumericallyIntegrable(x->exp(-x^4), (-Inf,Inf))

# ╔═╡ c2c72685-1d11-4691-bd00-5c9835fc88cb
g2 = NumericallyIntegrable(x->exp(-x^4), (-5,5))

# ╔═╡ d3c70173-f5c4-47e4-805c-f1cae4832014
begin
	plot()
	stephist!(rand(g, 100_000))
	stephist!(rand(g2, 100_000))
end

# ╔═╡ fff1bb78-e681-4f41-982f-48e0505d2c87
x_range = (2.0, 6.0)

# ╔═╡ 6d3f5c0b-43a8-42f6-aa7f-7a39d4a8f9b0
begin
	model_bkg(x) = truncated(Exponential(x.τ), x_range...)
	model_sig(x) = Normal(x.m, x.σ)
	function model_tot(x)
		fb = 1-x.f1-x.f2-x.f3-x.f4
		MixtureModel(
		[model_sig(x.s1),
			model_sig(x.s2),
			model_sig(x.s3),
			model_sig(x.s4),
			model_bkg(x.bkg)], [x.f1, x.f2, x.f3, x.f4, fb])
	end
end

# ╔═╡ 7587d930-433f-4c8b-81c3-2d8a6c59f6d0
pars = ComponentVector(
    s1=ComponentVector(m=3.2, σ=0.1),
    s2=ComponentVector(m=3.4, σ=0.04),
    s3=ComponentVector(m=4.6, σ=0.03),
    s4=ComponentVector(m=4.9, σ=0.051),
	# 
    bkg=ComponentVector(τ=2),
    f1=0.02,
    f2=0.03,
    f3=0.04,
    f4=0.03
)

# ╔═╡ a52ca9c2-c413-430e-9951-10fe8792b245
model = model_tot(pars)

# ╔═╡ ed4e775b-e0f9-47a6-8262-ab15a82f1143
let
	bins=range(x_range..., 30)
	dx = bins[2]-bins[1]
	plot()
	t = model_tot(pars)
	plot!(x->pdf(t,x), x_range..., lw=2, lc=:red, lab="Total")
	# 
	b = t.components[end]
	fb = t.prior.p[end]
	plot!(x->pdf(b,x) * fb, x_range..., lab="Background")
	#
	for i in 1:4
		s = t.components[i]
		fs = t.prior.p[i]
		plot!(x->pdf(s,x) * fs, x_range..., fill=0, alpha=0.4, lab="s$i")
	end
	plot!()
end

# ╔═╡ 7aaa9bba-a41c-4d98-97cc-ad1aac6fe078
data = rand(model, 10_000)

# ╔═╡ 1437790b-d282-4e72-b0a7-c789383a2cda
h = Hist1D(data, binedges=range(2,6,100))

# ╔═╡ ee58843a-bf95-4f0e-8da6-21584d11c9d6
const nSamples = 1_000

# ╔═╡ a4b4df4e-1f21-4542-a54a-2fe9ad3ac91a
const m, Γ = 0.77, 0.15

# ╔═╡ c29bf811-7e3a-414c-a2a4-b8d9af80482f
let
	scatter(bincenters(h), bincounts(h), yerr=sqrt.(bincounts(h)))
	f(x) = pdf(model, x) * integral(h; width=true)
	plot!(f, h.binedges[1][1], h.binedges[1][end])
end

# ╔═╡ f1f7e74d-f4e1-4fd4-a28b-a99d35cf294d
function chi2(h, d::UnivariateDistribution)
	scale = integral(h; width=true)
	f(x) = pdf(d, x) * scale
	chi2(h, f)
end

# ╔═╡ 92eccb71-70cd-48c3-bda3-188629ac6005
function chi2(h, model)
	# 
	xv = bincenters(h)
	yv = bincounts(h)
	δyv = sqrt.(yv)
	# 
	yv_pred = model.(xv)
	chi2 = sum(@. (yv_pred - yv)^2 / yv)
	# 
	return chi2
end

# ╔═╡ 3363d8d2-7dab-4034-ba23-8c77cfa21131
chi2_ndf = map(1:nSamples) do _
	data = rand(model, 10_000)
	h = Hist1D(data, binedges=range(2,6,100))
	chi2(h, model) / nbins(h)
end

# ╔═╡ 92f2e657-b389-4743-ab63-a71d6acaf39f
stephist(chi2_ndf, bins=100)

# ╔═╡ Cell order:
# ╠═a4f34cb4-170b-11f0-1df4-114a305c6853
# ╠═647abad9-5891-4154-8f91-2163b5825434
# ╠═29e9dda3-6b55-4dc7-af7e-f2d0b20bd3ae
# ╠═c2c72685-1d11-4691-bd00-5c9835fc88cb
# ╠═d3c70173-f5c4-47e4-805c-f1cae4832014
# ╠═fff1bb78-e681-4f41-982f-48e0505d2c87
# ╠═6d3f5c0b-43a8-42f6-aa7f-7a39d4a8f9b0
# ╠═7587d930-433f-4c8b-81c3-2d8a6c59f6d0
# ╠═a52ca9c2-c413-430e-9951-10fe8792b245
# ╠═ed4e775b-e0f9-47a6-8262-ab15a82f1143
# ╠═7aaa9bba-a41c-4d98-97cc-ad1aac6fe078
# ╠═1437790b-d282-4e72-b0a7-c789383a2cda
# ╠═ee58843a-bf95-4f0e-8da6-21584d11c9d6
# ╠═3363d8d2-7dab-4034-ba23-8c77cfa21131
# ╠═a4b4df4e-1f21-4542-a54a-2fe9ad3ac91a
# ╠═92f2e657-b389-4743-ab63-a71d6acaf39f
# ╠═c29bf811-7e3a-414c-a2a4-b8d9af80482f
# ╠═f1f7e74d-f4e1-4fd4-a28b-a99d35cf294d
# ╠═92eccb71-70cd-48c3-bda3-188629ac6005
