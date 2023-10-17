### A Pluto.jl notebook ###
# v0.19.29

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
        Pkg.PackageSpec("NLsolve"),
		Pkg.PackageSpec("ForwardDiff")
    ])
    # 
    using RecipesBase
    using ThreeBodyDecay
    using Plots
    using Polynomials
    using Parameters
    using NLsolve
	using ForwardDiff
	import ForwardDiff: Dual, partials, value, gradient, jacobian
	using LinearAlgebra
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

# ╔═╡ 6ce6fa9b-a227-40e4-b12a-6db062078fa1
begin
	import ThreeBodyDecay: Invariants
	# 
	function Invariants(ms::NamedTuple{(:m1, :m2, :m3, :m0)}; σ1=-1.0, σ2=-1.0, σ3=-1.0)
		sign(σ1) + sign(σ2) + sign(σ3) != 1 && error("the method works with TWO invariants given: $((σ1,σ2,σ3))")
		
		ms2 = ms |> collect .|> abs2
	    σ3 < 0 && return NamedTuple{(:σ1, :σ2, :σ3)}((σ1, σ2, σ3=sum(ms2) - σ1 - σ2))
	    σ1 < 0 && return NamedTuple{(:σ1, :σ2, :σ3)}((sum(ms2) - σ2 - σ3, σ2, σ3))
	    return NamedTuple{(:σ1, :σ2, :σ3)}((σ1=σ1, σ2=sum(ms2) - σ3 - σ1, σ3=σ3))
	end
	# 
	function _Invariants(ms::NamedTuple{(:m1, :m2, :m3, :m0)}; σ1=-1.0, σ2=-1.0, σ3=-1.0)
		sign(σ1) + sign(σ2) + sign(σ3) != 1 && error("the method works with TWO invariants given: $((σ1,σ2,σ3))")
		
		ms2 = ms |> collect .|> abs2
	    σ3 < 0 && return NamedTuple{(:σ1, :σ2, :σ3)}((σ1, σ2, σ3=sum(ms2) - σ1 - σ2))
	    σ1 < 0 && return NamedTuple{(:σ1, :σ2, :σ3)}((sum(ms2) - σ2 - σ3, σ2, σ3))
	    return NamedTuple{(:σ1, :σ2, :σ3)}((σ1=σ1, σ2=sum(ms2) - σ3 - σ1, σ3=σ3))
	end
end

# ╔═╡ 823b5793-8c1a-4277-95ed-0f84faf5f088
function alwaysin(ms)
    @unpack m0, m1, m2, m3 = ms
    # 
	ms2 = ms |> collect .|> abs2
    σ1 = ((m2 + m3 + m0 - m1) / 2)^2
    σ3 = σ3of1(0.1, σ1, ms2)
    σs2 = _Invariants(ms; σ1, σ3)
    # 
    σ2 = ((m3 + m1 + m0 - m2) / 2)^2
    σ1 = σ1of2(0, σ2, ms2)
    σs3 = _Invariants(ms; σ2, σ1)
    # 
    σ3 = ((m1 + m2 + m0 - m3) / 2)^2
    σ2 = σ2of3(0, σ3, ms2)
    σs1 = _Invariants(ms; σ3, σ2)

    σs = NamedTuple{(:σ1, :σ2, :σ3)}(sum((σs1, σs2, σs3) .|> collect) / 3)
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
begin
	σs_P(θ, σs0) =
	    Polynomial([σs0[1], -cos(θ)]),
	    Polynomial([σs0[2], cos(θ + π / 3)]),
	    Polynomial([σs0[3], cos(θ - π / 3)])
	
	function rborder(θ, ms)
		σs0 = alwaysin(ms)[end]
	    # 
	    ϕ = Base.Fix2(Kibble, ms^2)
		return minimum(filter(x -> x > 0, roots(ϕ(σs_P(θ, σs0)))))
	end
	function σsborder(θ, ms)
		σs0 = alwaysin(ms)[end]
		# 
	    map(σs_P(θ, σs0)) do P
			P(rborder(θ, ms))
	    end
	end
	
	function fixedborder(ms::ThreeBodyDecay.MassTuple; Nx::Int=300)
	    θs = range(-π, π, length=Nx)
	    return ThreeBodyDecay.MandestamTuple.(σsborder.(θs, Ref(ms)))
	end
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

# ╔═╡ 7f6e8261-aba4-4d62-b6e9-4454e5a00a14
md"""
## The $(\alpha,\theta)$ mapping
"""

# ╔═╡ 8cf51a81-aae5-4680-b900-968d0b3df5ff
function scalledmasses(α, ms)
    m_min = minimum([ms.m1, ms.m2, ms.m3])
    return (;
        m1=m_min + α * (ms.m1 - m_min),
        m2=m_min + α * (ms.m2 - m_min),
        m3=m_min + α * (ms.m3 - m_min),
        m0=3m_min + α * (ms.m0 - 3m_min))
end

# ╔═╡ 44efb115-adde-4c67-9833-6451b9952220
rθborder(ms) = σs2rθ.(fixedborder(ms), Ref(ms))

# ╔═╡ 5c0adba5-6cb3-4fe6-8f44-ebeb9cca26b4
function rθ2σs(θr::NamedTuple{(:θ, :r)},
    ms::NamedTuple{(:m1, :m2, :m3, :m0)})
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
    NamedTuple{(:σ1, :σ2, :σ3)}(σs)
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

# ╔═╡ fcf282e9-3c30-4ccb-99e0-f8c8d5ef0d8c
md"""
## Transformation to alpha
"""

# ╔═╡ 0124a4d4-4aa3-4875-b47f-eb3c2cada230
function Kibble_α(sqrtα, x_rθ, ms)
    msα = scalledmasses(sqrtα^2, ms)
	rθ = NamedTuple{(:θ,:r)}(x_rθ)
    Kibble(rθ2σs(rθ, msα), msα |> collect .|> abs2)
end

# ╔═╡ 33eb8faa-e1b0-4126-bc47-6e2aee77c44d
rθ2sqrtα(x_rθ, ms) = nlsolve(x->Kibble_α(x[1], x_rθ, ms), [1.0]).zero[1]

# ╔═╡ 86015c97-0ecd-4196-9636-35ad4129c753
function rθ2sqrtα(d::Vector{Dual{T,V,N}}, ms) where {T,V,N}
	x0 = value.(d)
	y0 = rθ2sqrtα(x0, ms)
	# 
	dfdy = ForwardDiff.derivative(sqrtα -> Kibble_α(sqrtα, x0, ms), y0)
    dfdx = ForwardDiff.gradient(x->Kibble_α(y0, x, ms), x0)
    dydx = -inv(dfdy) * dfdx
	# 
	nt = sum(i->dydx[i]*partials.(d, i), 1:N)
	Dual{T}(y0, Tuple(nt))
end

# ╔═╡ aa661d42-5273-4f63-a71a-795a20d1da14
function σs2ζχ(σs, ms)
    rθ = σs2rθ(σs, ms)
	x_rθ = rθ |> collect
	# 
	sqrtα = rθ2sqrtα(x_rθ, ms)
    α = sqrtα^2
    α > 1 && error("α = $α > 1")
	# 
	θ = rθ.θ
    [α*cos(θ), α*sin(θ)]
end

# ╔═╡ af7c35b6-971f-4860-811b-4a6b9acb4017
function ζχ2σs(ζχ, ms)
	ζ,χ = ζχ
	α = norm(ζχ)
	α > 1 && error("α = $α > 1") 
	θ = atan(χ,ζ)
	# 
	σs0 = alwaysin(ms)[end]
	r = rborder(θ, scalledmasses(α, ms))
	# 
	σsv = map(σs_P(θ, σs0)) do P
		P(r)
	end
	return NamedTuple{(:σ1, :σ2, :σ3)}(σsv)
end

# ╔═╡ aab3ddac-fae2-4a68-afe4-b8828df3c237
ζχ2σs(σs2ζχ(σs_i, ms), ms), σs_i

# ╔═╡ 49ae880a-8acb-425c-9031-ad94a952ef30
xxx(σ1σ2) = σs2ζχ(_Invariants(ms; σ1=σ1σ2[1], σ2 = σ1σ2[2]), ms)

# ╔═╡ 25c6c066-e511-459a-990f-7d5063f23483
σs2ζχ_jacobian(σs, ms) = ForwardDiff.jacobian(xxx, collect(σs)[1:2])

# ╔═╡ 2257dd4f-9480-404b-a43f-6fdb3ab2a130
md"""
## Cross check
"""

# ╔═╡ 58d6ef6b-855c-4057-916d-97033c102b2f
jaco = σs2ζχ_jacobian([110,160], ms)

# ╔═╡ 14c7cdd4-aaef-4259-a263-7ef9f309e4cd
grads = hcat(
	ForwardDiff.gradient([110, 160]) do σ1σ2
			xxx(σ1σ2)[1]
	end,
	ForwardDiff.gradient([110, 160]) do σ1σ2
		xxx(σ1σ2)[2]
	end) |> transpose

# ╔═╡ 5c03220f-01d7-4be5-8f09-97dff3f5a918
ForwardDiff.gradient(σ1σ2->xxx(σ1σ2)[1], [110, 160])

# ╔═╡ 99a405a0-9062-4f69-b9fa-7428c07279d1
let ϵ = 1e-6
	(xxx([110,160+ϵ])[1] - xxx([110,160])[1])/ϵ
end

# ╔═╡ 7e174855-5875-48ea-b51a-227d13b9fb5f
finite = let ϵ = 1e-6
	hcat(
		(xxx([110+ϵ,160]) - xxx([110-ϵ,160])) / 2ϵ,
		(xxx([110,160+ϵ]) - xxx([110,160-ϵ])) / 2ϵ
	) |> transpose
end

# ╔═╡ e01118fb-b988-43ea-8fc7-973d2dd1789d
yyy(ζχ) = let ϵ = 1e-6
	σ12 = (ζχ2σs(ζχ, ms) |> collect)[1:2]
	hcat(
		(xxx(σ12 + [0,ϵ]) - xxx(σ12)) / 2ϵ,
		(xxx(σ12 + [ϵ,0]) - xxx(σ12)) / 2ϵ
	) |> det
end

# ╔═╡ fdb11760-24c1-475d-b1d6-373ff907a9db
yyy([0.3,0.5])

# ╔═╡ 31903638-516e-4562-9515-1c371810149c
finite ./ grads

# ╔═╡ 6300912a-3e8e-4868-9944-1e7bd64cb5b9
ζχ_jacobian(ζχ, ms) = det(σs2ζχ_jacobian(ζχ2σs(ζχ, ms), ms))

# ╔═╡ 22378478-cf42-49c0-aa9b-557e3f036d7e
ζχ_jacobian([0.2,0.3], ms)

# ╔═╡ a2efb67a-2551-490d-b36b-dc2fb1816bdb
let
	xv = range(-0.1, 0.1, 50)
	yv = range(-0.1, 0.1, 50) # ζχ_jacobian([x,y], ms)
	calv = [norm([x,y]) > 1 ? NaN : yyy([x,y]) for x in xv, y in yv]
	heatmap(xv,yv,calv, aspect_ratio=1, c=:viridis, colorbar=true)
end

# ╔═╡ 55214e81-d896-4b58-b2c1-75c84ea35ffd
begin
	function σs2α(rθ, ms)
	    sqrtα = nlsolve(x->Kibble_α(x[1], rθ, ms), [1.0]).zero[1]
	    α = sqrtα^2
	    α > 1 && error("α = $α > 1")
	    return sqrtα^2
	end
	function σs2θα(σs, ms)
	    rθ = σs2rθ(σs, ms)
	    α = σs2α(rθ, ms)
	    (; rθ.θ, α=α)
	end
end

# ╔═╡ 1dcc0b4b-cff0-4448-bf89-e32dcd4ceda8
data_σs = flatDalitzPlotSample(ms; Nev=10_000);

# ╔═╡ d0dc58cc-d818-4e8e-b961-6d39bea67691
data_ζχ = Tuple.(σs2ζχ.(data_σs, Ref(ms)))

# ╔═╡ d0e90d94-f209-4473-93b3-2fc20ead1e68
histogram2d(data_ζχ, aspect_ratio=1, bins=100)

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
# ╠═6ce6fa9b-a227-40e4-b12a-6db062078fa1
# ╠═823b5793-8c1a-4277-95ed-0f84faf5f088
# ╟─7f6e8261-aba4-4d62-b6e9-4454e5a00a14
# ╠═8cf51a81-aae5-4680-b900-968d0b3df5ff
# ╠═44efb115-adde-4c67-9833-6451b9952220
# ╠═1bab77b3-9efb-425c-9227-6b6312f3764f
# ╠═5c0adba5-6cb3-4fe6-8f44-ebeb9cca26b4
# ╟─0f94bc6d-0084-4136-8613-de48bc877149
# ╠═89b5cf5e-20bf-4de0-9522-f019536e84b6
# ╟─fcf282e9-3c30-4ccb-99e0-f8c8d5ef0d8c
# ╠═0124a4d4-4aa3-4875-b47f-eb3c2cada230
# ╠═33eb8faa-e1b0-4126-bc47-6e2aee77c44d
# ╠═86015c97-0ecd-4196-9636-35ad4129c753
# ╠═aa661d42-5273-4f63-a71a-795a20d1da14
# ╠═af7c35b6-971f-4860-811b-4a6b9acb4017
# ╠═aab3ddac-fae2-4a68-afe4-b8828df3c237
# ╠═49ae880a-8acb-425c-9031-ad94a952ef30
# ╠═25c6c066-e511-459a-990f-7d5063f23483
# ╟─2257dd4f-9480-404b-a43f-6fdb3ab2a130
# ╠═58d6ef6b-855c-4057-916d-97033c102b2f
# ╠═14c7cdd4-aaef-4259-a263-7ef9f309e4cd
# ╠═5c03220f-01d7-4be5-8f09-97dff3f5a918
# ╠═99a405a0-9062-4f69-b9fa-7428c07279d1
# ╠═7e174855-5875-48ea-b51a-227d13b9fb5f
# ╠═e01118fb-b988-43ea-8fc7-973d2dd1789d
# ╠═fdb11760-24c1-475d-b1d6-373ff907a9db
# ╠═31903638-516e-4562-9515-1c371810149c
# ╠═6300912a-3e8e-4868-9944-1e7bd64cb5b9
# ╠═22378478-cf42-49c0-aa9b-557e3f036d7e
# ╠═a2efb67a-2551-490d-b36b-dc2fb1816bdb
# ╠═55214e81-d896-4b58-b2c1-75c84ea35ffd
# ╠═1dcc0b4b-cff0-4448-bf89-e32dcd4ceda8
# ╠═d0dc58cc-d818-4e8e-b961-6d39bea67691
# ╠═d0e90d94-f209-4473-93b3-2fc20ead1e68
