### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ b2323759-cbb9-4a72-ac40-46076819274f
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"),
		Pkg.PackageSpec("Polynomials"),
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("Parameters"),
		Pkg.PackageSpec("RecipesBase"),
		Pkg.PackageSpec("NLsolve")
		])
	# 
	using RecipesBase
	using ThreeBodyDecay
	using Plots
	using Polynomials
	using Parameters
	using NLsolve
end

# ╔═╡ edd93652-4d2d-4b50-a8f2-61eea1ea9695
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lw=1.2, lab="", colorbar=false)

# ╔═╡ 99a30a66-2ae7-4914-9e68-b97f95dea378
function rθdalitz(σs, ms)
    thresholds = ((ms.m2 + ms.m3), (ms.m3 + ms.m1), (ms.m1 + ms.m2))
    T0 = sum(abs2, ms) .- sum(abs2, thresholds)
    Δσs = Tuple(σs) .- thresholds .^ 2
    # 
    rcosθ = T0 / 3 - Δσs[1]
    rcosθmπ3 = Δσs[3] - T0 / 3
    rsinθ = 2 / √3 * (rcosθmπ3 - rcosθ / 2)
    # 
    θ = atan(rsinθ, rcosθ)
    r = rcosθ / cos(θ)
    # 
    return (; θ, r)
end

# ╔═╡ f59fbbff-8470-4926-89c3-14dd5ab33ee8
ms = ThreeBodyMasses(2, 2, 6; m0=14)

# ╔═╡ 4bf79a8c-8f2e-4561-b886-1964ba95a11b
σsv = border(ms)

# ╔═╡ ce4ca0e9-ee3f-47fc-a56c-570c33d61a9a
rθv = rθdalitz.(σsv, Ref(ms))

# ╔═╡ 947586cf-63a9-4d5e-9caf-a175e3ec8bb2
plot(
	plot(NamedTuple{(:σ1, :σ2)}.(σsv), title="Cartesian Dalitz plot"),
	plot(rθv, title="Polar dalitz plot", proj=:polar, lims=(0, :auto))
)

# ╔═╡ ab731426-ce4b-4e8f-bb67-a7b9391269e0
function rlimsdalitz(θ, ms)
    thresholds = ((ms.m2 + ms.m3), (ms.m3 + ms.m1), (ms.m1 + ms.m2))
    T0 = sum(abs2, ms) .- sum(abs2, thresholds)
    # 
    ϕ(σs) = Kibble(σs, ms^2)
    σs = thresholds .^ 2 .+ ThreeBodyDecay.polardalitz2invariants(θ, T0)
    rmax = minimum(filter(x -> x > 0, roots(ϕ(σs))))
    return rmax
end

# ╔═╡ 3ff2fbbf-808f-4cde-ac3f-6e4ebf265df9
md"""
## $(r,\theta)$ mapping
"""

# ╔═╡ 8edcf60a-f68e-4419-a45b-717ca8a071ad
let
    θv = range(-π, π, length=200)
    rv = range(0, 1, length=101)
    zv = [rlimsdalitz(θ, ms) for r in rv, θ in θv]
    Plots.heatmap(θv, rv, zv, proj=:polar, title="Phase-space intensity")
end

# ╔═╡ fa2e688f-c434-4e1e-8b5e-c58cf6ccf30f
md"""
The point of $r=0$ is not holomorphic of the scalaring transformation
"""

# ╔═╡ 7f6e8261-aba4-4d62-b6e9-4454e5a00a14
md"""
## The $(\alpha,\theta)$ mapping
"""

# ╔═╡ 8cf51a81-aae5-4680-b900-968d0b3df5ff
function scalledmasses(α,ms)
	m_min = minimum([ms.m1,ms.m2,ms.m3])
	return ThreeBodyMasses(
		m_min+α*(ms.m1-m_min),
		m_min+α*(ms.m2-m_min),
		m_min+α*(ms.m3-m_min);
		m0=3m_min+α*(ms.m0-3m_min))
end

# ╔═╡ 44efb115-adde-4c67-9833-6451b9952220
rθborder(ms) = rθdalitz.(border(ms), Ref(ms))

# ╔═╡ 5c0adba5-6cb3-4fe6-8f44-ebeb9cca26b4
function rθ2σs(θr::NamedTuple{(:θ,:r)},
	ms::ThreeBodyDecay.MassTuple)
	#
	@unpack θ,r = θr
	# 
	thresholds = ((ms.m2 + ms.m3), (ms.m3 + ms.m1), (ms.m1 + ms.m2))
    T0 = sum(abs2, ms) .- sum(abs2, thresholds)
	# 
	Ps = ThreeBodyDecay.polardalitz2invariants(θ, T0)
	σs = thresholds .^ 2 .+ map.(Ps, r)
	# 
	ThreeBodyDecay.MandestamTuple(σs)
end

# ╔═╡ 0f94bc6d-0084-4136-8613-de48bc877149
md"""
Test that the transformation ∘ inv_transformation = 1
"""

# ╔═╡ 89b5cf5e-20bf-4de0-9522-f019536e84b6
begin
	σs_i = randomPoint(ms)
	σs_i_back = rθ2σs(rθdalitz(σs_i,ms), ms)
	@assert prod(collect(σs_i_back) .≈ collect(σs_i))
end

# ╔═╡ 1bab77b3-9efb-425c-9227-6b6312f3764f
begin
	plot(rθborder(ms), proj =:polar, lims=(0, :auto))
	for α  in 0.1:0.1:0.9
		plot!(rθborder(scalledmasses(α,ms)), proj =:polar, lims=(0, :auto))
	end
	scatter!(rθdalitz.([σs_i], Ref(ms)), proj =:polar, m=(6,:red))
	plot!()
end

# ╔═╡ 55214e81-d896-4b58-b2c1-75c84ea35ffd
function σs2α(rθ, ms)
	function Kibble_α(sqrtα)
		msα = scalledmasses(sqrtα^2,ms)
		Kibble(rθ2σs(rθ, msα), msα^2)
	end
	#
	α0 = rθ.r / rlimsdalitz(rθ.θ, ms)
	sqrtα = nlsolve(n_ary(Kibble_α), [sqrt(α0)]).zero[1]
	α = sqrtα^2
	α > 1 && error("α = $α > 1")
	return sqrtα^2
end

# ╔═╡ 24fee2f9-c808-41e3-ba6c-aa6824fd3174
function σs2θα(σs, ms)
	rθ = rθdalitz(σs, ms)
	α = σs2α(rθ, ms)
	(; rθ.θ, α=α)
end

# ╔═╡ 1dcc0b4b-cff0-4448-bf89-e32dcd4ceda8
data_σs = flatDalitzPlotSample(ms; Nev=100_000);

# ╔═╡ efc7c156-716c-4725-b62b-76682bdd3498
data_σs_selected = filter(data_σs) do σs
	# rθdalitz(σs, ms).r < 10
	true
end;

# ╔═╡ 5d88295e-5209-452c-9c9e-d2cf9eeabd32
data_θα = σs2θα.(data_σs_selected, Ref(ms));

# ╔═╡ 48353bd2-ddb7-4eb1-858d-ed6406b93da8
data_xy = map(data_θα) do (θ,r)
	(r*cos(θ), r*sin(θ))
end;

# ╔═╡ 736d2191-3682-4ca5-989d-099c13393132
histogram2d(data_xy, aspect_ratio=1, bins=100)

# ╔═╡ Cell order:
# ╠═b2323759-cbb9-4a72-ac40-46076819274f
# ╠═edd93652-4d2d-4b50-a8f2-61eea1ea9695
# ╠═99a30a66-2ae7-4914-9e68-b97f95dea378
# ╠═f59fbbff-8470-4926-89c3-14dd5ab33ee8
# ╠═4bf79a8c-8f2e-4561-b886-1964ba95a11b
# ╠═ce4ca0e9-ee3f-47fc-a56c-570c33d61a9a
# ╠═947586cf-63a9-4d5e-9caf-a175e3ec8bb2
# ╠═ab731426-ce4b-4e8f-bb67-a7b9391269e0
# ╟─3ff2fbbf-808f-4cde-ac3f-6e4ebf265df9
# ╠═8edcf60a-f68e-4419-a45b-717ca8a071ad
# ╟─fa2e688f-c434-4e1e-8b5e-c58cf6ccf30f
# ╟─7f6e8261-aba4-4d62-b6e9-4454e5a00a14
# ╠═8cf51a81-aae5-4680-b900-968d0b3df5ff
# ╠═44efb115-adde-4c67-9833-6451b9952220
# ╠═1bab77b3-9efb-425c-9227-6b6312f3764f
# ╠═5c0adba5-6cb3-4fe6-8f44-ebeb9cca26b4
# ╟─0f94bc6d-0084-4136-8613-de48bc877149
# ╠═89b5cf5e-20bf-4de0-9522-f019536e84b6
# ╠═55214e81-d896-4b58-b2c1-75c84ea35ffd
# ╠═24fee2f9-c808-41e3-ba6c-aa6824fd3174
# ╠═1dcc0b4b-cff0-4448-bf89-e32dcd4ceda8
# ╠═efc7c156-716c-4725-b62b-76682bdd3498
# ╠═5d88295e-5209-452c-9c9e-d2cf9eeabd32
# ╠═48353bd2-ddb7-4eb1-858d-ed6406b93da8
# ╠═736d2191-3682-4ca5-989d-099c13393132
