### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ b2323759-cbb9-4a72-ac40-46076819274f
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"),
		Pkg.PackageSpec("Polynomials"),
		Pkg.PackageSpec("Plots")
		])
	# 
	using ThreeBodyDecay
	using Plots
	using Polynomials
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
    return (; r, θ)
end

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

# ╔═╡ f59fbbff-8470-4926-89c3-14dd5ab33ee8
ms = ThreeBodyMasses(1, 2, 4; m0=13)

# ╔═╡ 4bf79a8c-8f2e-4561-b886-1964ba95a11b
σsv = border(ms)

# ╔═╡ ce4ca0e9-ee3f-47fc-a56c-570c33d61a9a
rθv = rθdalitz.(σsv, Ref(ms))

# ╔═╡ 947586cf-63a9-4d5e-9caf-a175e3ec8bb2
plot(
	plot(NamedTuple{(:σ1, :σ2)}.(σsv), title="Cartesian Dalitz plot"),
	plot(title="Polar dalitz plot",
	    getproperty.(rθv, :θ),
	    getproperty.(rθv, :r), proj=:polar, lims=(0, :auto))
)


# ╔═╡ 8edcf60a-f68e-4419-a45b-717ca8a071ad
let
    θv = range(-π, π, length=200)
    rv = range(0, 1, length=101)
    zv = [rlimsdalitz(θ, ms) for r in rv, θ in θv]
    Plots.heatmap(θv, rv, zv, proj=:polar, title="Phase-space intensity")
end

# ╔═╡ 46f70e9a-1818-4d6c-9bbc-fe1f23960482
begin
	dv = flatDalitzPlotSample(ms; Nev=50_000)
	rθdv = rθdalitz.(dv, Ref(ms))
end ;

# ╔═╡ 3d055f81-5fbb-4575-aaaa-8ebb0f22666c
scatter(
    getproperty.(rθdv, :θ),
    getproperty.(rθdv, :r) ./ rlimsdalitz.(getproperty.(rθdv, :θ), Ref(ms)), 
	proj=:polar, lims=(0, :auto), c=2, ms=0.1,
	title="Phase-space sample")

# ╔═╡ Cell order:
# ╠═b2323759-cbb9-4a72-ac40-46076819274f
# ╠═edd93652-4d2d-4b50-a8f2-61eea1ea9695
# ╠═99a30a66-2ae7-4914-9e68-b97f95dea378
# ╠═ab731426-ce4b-4e8f-bb67-a7b9391269e0
# ╠═f59fbbff-8470-4926-89c3-14dd5ab33ee8
# ╠═4bf79a8c-8f2e-4561-b886-1964ba95a11b
# ╠═ce4ca0e9-ee3f-47fc-a56c-570c33d61a9a
# ╠═947586cf-63a9-4d5e-9caf-a175e3ec8bb2
# ╠═8edcf60a-f68e-4419-a45b-717ca8a071ad
# ╠═46f70e9a-1818-4d6c-9bbc-fe1f23960482
# ╠═3d055f81-5fbb-4575-aaaa-8ebb0f22666c
