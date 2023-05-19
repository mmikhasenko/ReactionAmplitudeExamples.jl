### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ b2323759-cbb9-4a72-ac40-46076819274f
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(mktempdir())
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

# ╔═╡ 0270aecb-b6c2-4d29-8760-77642194838e
md"""
# Circular Dalitz Mapping

The notebook demonstrate the mapping algorithm of an arbitrary dalitz plot to a unite circle that preserve the symmatry between the Mandelsam variables
"""

# ╔═╡ edd93652-4d2d-4b50-a8f2-61eea1ea9695
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lw=1.2, lab="", colorbar=false)

# ╔═╡ f59fbbff-8470-4926-89c3-14dd5ab33ee8
ms = ThreeBodyMasses(3, 1, 9; m0=14)

# ╔═╡ 823b5793-8c1a-4277-95ed-0f84faf5f088
function alwaysin(ms)
    @unpack m0, m1, m2, m3 = ms
    # 
    σ1 = ((m2 + m3 + m0 - m1) / 2)^2
    σ3 = σ3of1(0.1, σ1, ms^2)
    σs2 = Invariants(ms; σ1, σ3)
    # 
    σ2 = ((m3 + m1 + m0 - m2) / 2)^2
    σ1 = σ1of2(0, σ2, ms^2)
    σs3 = Invariants(ms; σ2, σ1)
    # 
    σ3 = ((m1 + m2 + m0 - m3) / 2)^2
    σ2 = σ2of3(0, σ3, ms^2)
    σs1 = Invariants(ms; σ3, σ2)

    σs = ThreeBodyDecay.MandestamTuple(sum((σs1, σs2, σs3) .|> collect) / 3)
    [σs1, σs2, σs3, σs]
end

# ╔═╡ fc40eced-a66e-4e08-ad1c-9a5dc39d4ef5
function σs2rθ(σs, ms)
    #
    σs0 = alwaysin(ms)[end]
    Δσs0 = Tuple(σs0)
    Δσs = Tuple(σs)
    # 
    rcosθ = Δσs0[1] - Δσs[1]
    rcosθmπ3 = Δσs[3] - Δσs0[3]
    rsinθ = 2 / √3 * (rcosθmπ3 - rcosθ / 2)
    # 
    θ = atan(rsinθ, rcosθ)
    r = rcosθ / cos(θ)
    # 
    return (; θ, r)
end

# ╔═╡ dc1b8aa0-1cd6-44ce-9e3f-0d56bfc269ac
function fixedborder(ms::ThreeBodyDecay.MassTuple; Nx::Int=300)
    # 
    σs0 = alwaysin(ms)[end]
    # 
    σs_P(θ) =
        Polynomial([σs0[1], -cos(θ)]),
        Polynomial([σs0[2], cos(θ + π / 3)]),
        Polynomial([σs0[3], cos(θ - π / 3)])
    # 
    ϕ = Base.Fix2(Kibble, ms^2)
    rborder(θ) = minimum(filter(x -> x > 0, roots(ϕ(σs_P(θ)))))
    σsborder(θ) =
        map(σs_P(θ)) do P
            P(rborder(θ))
        end
    θs = range(-π, π, length=Nx)
    return ThreeBodyDecay.MandestamTuple.(σsborder.(θs))
end

# ╔═╡ 4bf79a8c-8f2e-4561-b886-1964ba95a11b
σsv = fixedborder(ms)

# ╔═╡ ce4ca0e9-ee3f-47fc-a56c-570c33d61a9a
rθv = σs2rθ.(σsv, Ref(ms))

# ╔═╡ 947586cf-63a9-4d5e-9caf-a175e3ec8bb2
plot(
    plot(NamedTuple{(:σ1, :σ2)}.(σsv), title="Cartesian Dalitz plot"),
    plot(rθv, title="Polar dalitz plot", proj=:polar, lims=(0, :auto))
)

# ╔═╡ ab731426-ce4b-4e8f-bb67-a7b9391269e0
function rlimsdalitz(θ, ms)
    σs0 = alwaysin(ms)[end]
    # 
    σs_P =
        Polynomial([σs0[1], -cos(θ)]),
        Polynomial([σs0[2], cos(θ + π / 3)]),
        Polynomial([σs0[3], cos(θ - π / 3)])
    # 
    ϕ = Base.Fix2(Kibble, ms^2)
    rmax = minimum(filter(x -> x > 0, roots(ϕ(σs_P))))
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
function scalledmasses(α, ms)
    m_min = minimum([ms.m1, ms.m2, ms.m3])
    return ThreeBodyMasses(
        m_min + α * (ms.m1 - m_min),
        m_min + α * (ms.m2 - m_min),
        m_min + α * (ms.m3 - m_min);
        m0=3m_min + α * (ms.m0 - 3m_min))
end

# ╔═╡ 44efb115-adde-4c67-9833-6451b9952220
rθborder(ms) = σs2rθ.(fixedborder(ms), Ref(ms))

# ╔═╡ 5c0adba5-6cb3-4fe6-8f44-ebeb9cca26b4
function rθ2σs(θr::NamedTuple{(:θ, :r)},
    ms::ThreeBodyDecay.MassTuple)
    #
    @unpack θ, r = θr
    # 
    σs0 = alwaysin(ms)[end]
    # 
    Ps =
        Polynomial([σs0[1], -cos(θ)]),
        Polynomial([σs0[2], cos(θ + π / 3)]),
        Polynomial([σs0[3], cos(θ - π / 3)])
    # 
    σs = map.(Ps, r)
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
    σs_i_back = rθ2σs(σs2rθ(σs_i, ms), ms)
    @assert prod(collect(σs_i_back) .≈ collect(σs_i))
end

# ╔═╡ 1bab77b3-9efb-425c-9227-6b6312f3764f
begin
    plot(rθborder(ms), proj=:polar, lims=(0, :auto))
    for α in 0.1:0.1:0.9
        plot!(rθborder(scalledmasses(α, ms)), proj=:polar, lims=(0, :auto))
    end
    scatter!(σs2rθ.([σs_i], Ref(ms)), proj=:polar, m=(6, :red))
    plot!()
    # 
    scatter!(σs2rθ.(alwaysin(ms)[1:3], Ref(ms)), proj=:polar, α=0.2)
    scatter!(σs2rθ.([alwaysin(ms)[end]], Ref(ms)), proj=:polar)
end

# ╔═╡ 55214e81-d896-4b58-b2c1-75c84ea35ffd
function σs2α(rθ, ms)
    function Kibble_α(sqrtα)
        msα = scalledmasses(sqrtα^2, ms)
        Kibble(rθ2σs(rθ, msα), msα^2)
    end
    #
    sqrtα = nlsolve(n_ary(Kibble_α), [1.0]).zero[1]
    α = sqrtα^2
    α > 1 && error("α = $α > 1")
    return sqrtα^2
end

# ╔═╡ 24fee2f9-c808-41e3-ba6c-aa6824fd3174
function σs2θα(σs, ms)
    rθ = σs2rθ(σs, ms)
    α = σs2α(rθ, ms)
    (; rθ.θ, α=α)
end

# ╔═╡ 1dcc0b4b-cff0-4448-bf89-e32dcd4ceda8
data_σs = flatDalitzPlotSample(ms; Nev=100_000);

# ╔═╡ 5d88295e-5209-452c-9c9e-d2cf9eeabd32
data_θα = σs2θα.(data_σs, Ref(ms));

# ╔═╡ 48353bd2-ddb7-4eb1-858d-ed6406b93da8
data_xy = map(data_θα) do (θ, r)
    (r * cos(θ), r * sin(θ))
end;

# ╔═╡ 736d2191-3682-4ca5-989d-099c13393132
histogram2d(data_xy, aspect_ratio=1, bins=100)

# ╔═╡ Cell order:
# ╟─0270aecb-b6c2-4d29-8760-77642194838e
# ╠═b2323759-cbb9-4a72-ac40-46076819274f
# ╠═edd93652-4d2d-4b50-a8f2-61eea1ea9695
# ╠═fc40eced-a66e-4e08-ad1c-9a5dc39d4ef5
# ╠═dc1b8aa0-1cd6-44ce-9e3f-0d56bfc269ac
# ╠═f59fbbff-8470-4926-89c3-14dd5ab33ee8
# ╠═4bf79a8c-8f2e-4561-b886-1964ba95a11b
# ╠═ce4ca0e9-ee3f-47fc-a56c-570c33d61a9a
# ╠═947586cf-63a9-4d5e-9caf-a175e3ec8bb2
# ╠═ab731426-ce4b-4e8f-bb67-a7b9391269e0
# ╠═823b5793-8c1a-4277-95ed-0f84faf5f088
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
# ╠═5d88295e-5209-452c-9c9e-d2cf9eeabd32
# ╠═48353bd2-ddb7-4eb1-858d-ed6406b93da8
# ╠═736d2191-3682-4ca5-989d-099c13393132
