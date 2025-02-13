### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 1c6e06ae-5747-11ec-0f31-9347d77fe04b
# ╠═╡ show_logs = false
begin
    import Pkg
    Pkg.activate(mktempdir())
    # 
    Pkg.add("AlgebraPDF")
    Pkg.add("QuadGK")
    Pkg.add("Plots")
    Pkg.add("LaTeXStrings")
    Pkg.add("PlutoUI")
    # 
    using AlgebraPDF
    using Plots
    using QuadGK
    using LaTeXStrings
	using PlutoUI
end

# ╔═╡ ddfc8f61-e540-4db6-8ee6-e5cce754066f
theme(:boxed)

# ╔═╡ 18663f9a-5d04-4530-84e7-cadf474643da
const sth = 0.6^2

# ╔═╡ ab563532-193c-4cdf-a837-13124ce4da78
md"""
## Unitarized background
$U(s) = B(s) + T(s)\left(c+ \frac{1}{\pi}\int_{\text{th}}^{\infty} \frac{B(s')\rho(s')}{s'-s} \mathrm{d}s'\right)$

with

$\text{Im}(T^{-1}) = -i\rho$
"""

# ╔═╡ 2fdf557e-d2bb-4739-9004-ced593298f09
begin
    struct UnitarizedBackground{B,T,R} <: AbstractFunctionWithParameters
        B::B
        T::T
        ρ::R
        ϵ::Float64
    end
    UnitarizedBackground(B, T, ρ) = UnitarizedBackground(B, T, ρ, 1e-6)
    #
    AlgebraPDF.pars(ub::UnitarizedBackground, isfree::Bool) =
        pars(ub.B, isfree) + pars(ub.T, isfree) + pars(ub.ρ, isfree)
    # 
    AlgebraPDF.func(ub::UnitarizedBackground, x::NumberOrTuple; p=pars(ub)) =
    #
        func(ub.B, x; p) + func(ub.T, x; p) *
                           quadgk(x′ -> func(ub.B, x′; p) * ub.ρ(x′) / (x′ - x - 1im * ub.ϵ), 0.6^2, Inf)[1] / π
end

# ╔═╡ b1577a00-af60-47c5-936e-1bc1134a133a
begin
    @makefuntype Background(s; p) = s < sth ? 0.0 : (s - sth)^p.α * exp(-p.β * s)
    b = Background((α=0.3, β=0.2)) * (b=2.2,)
end

# ╔═╡ e78d6e79-2d1f-4471-8489-1d7d8ed33dd8
plot(e -> b(e^2), 0.6, 2.5)

# ╔═╡ fbaa22eb-041a-4525-bd11-0b6ca0f0dc4f
begin
    @makefuntype Signal(s; p) = p.m₀ * p.Γ₀ / (p.m₀^2 - s - 1im * p.m₀ * p.Γ₀)
    sig = Signal((m₀=1.8, Γ₀=0.08))
end

# ╔═╡ 2041e2f1-87b7-4cbe-accd-c462a71212a9
plot(e -> abs2(sig(e^2)), 0.6, 2.5)

# ╔═╡ a4d25e58-2815-4610-b90f-2225023a4609
ρ(s) = 1.0

# ╔═╡ bb05f42f-d70b-4824-978a-8e398d164d3d
b_ub = UnitarizedBackground(b, sig, ρ)

# ╔═╡ 00374615-e9e5-4774-a8d9-2377416ce042
let
    plot()
    plot!(e -> abs2(b_ub(e^2)), 0.6, 2.5, lab="B+U(B)", fill=0, c=2, alpha=0.4, linealpha=1)
    plot!(e -> abs2(b(e^2)), 0.6, 2.5, lab="B", c=1, lw=2)
    # plot!(e -> abs2(sig(e^2)) * 70, 0.6, 2.5, lab="Signal")
end

# ╔═╡ d5182be7-b732-4b19-95d0-de1e4bc8160c
sub = sig * (c1=1.0+0.0im,) + b_ub * (c2=1.0,)

# ╔═╡ f25e816f-6e53-4250-a751-604f8fa463b4
@bind ϕ Slider(-π:0.01:π, default=0.0, show_value=true)

# ╔═╡ 6133607f-6972-41df-8f98-9c26694f7d4a
let
	c1 = 8.0*cis(ϕ)
	d = updatepar(sub, :c1, c1)
	# 
	plot(ylim=(0, 80))
	plot!(e -> abs2(b_ub(e^2)), 0.6, 2.5, lab="B+U(B)", c=2, fill=0, linealpha=1, alpha=0.4)
	plot!(e -> abs2(d(e^2)), 0.6, 2.5, lab="total", lw=2, c=1)
end

# ╔═╡ 5d08c920-c03d-4c43-a241-428856fa7faa
# ╠═╡ disabled = true
#=╠═╡
let  
	anim = @animate for ϕ in -π:0.3:π
		c1 = 8*cis(ϕ)
	    d = updatepar(sub, :c1, c1)
		# 
		plot(ylim=(0, 80))
	    plot!(e -> abs2(b_ub(e^2)), 0.6, 2.5, lab="B+U(B)", c=2, fill=0, linealpha=1, alpha=0.4)
	    plot!(e -> abs2(d(e^2)), 0.6, 2.5, lab="total", lw=2, c=1)
		# 
		vline!([sig.p.:m₀])
	end
	gif(anim, "test.gif", fps = 15)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═1c6e06ae-5747-11ec-0f31-9347d77fe04b
# ╠═ddfc8f61-e540-4db6-8ee6-e5cce754066f
# ╠═18663f9a-5d04-4530-84e7-cadf474643da
# ╠═b1577a00-af60-47c5-936e-1bc1134a133a
# ╠═e78d6e79-2d1f-4471-8489-1d7d8ed33dd8
# ╠═fbaa22eb-041a-4525-bd11-0b6ca0f0dc4f
# ╠═2041e2f1-87b7-4cbe-accd-c462a71212a9
# ╟─ab563532-193c-4cdf-a837-13124ce4da78
# ╠═2fdf557e-d2bb-4739-9004-ced593298f09
# ╠═a4d25e58-2815-4610-b90f-2225023a4609
# ╠═bb05f42f-d70b-4824-978a-8e398d164d3d
# ╠═00374615-e9e5-4774-a8d9-2377416ce042
# ╠═d5182be7-b732-4b19-95d0-de1e4bc8160c
# ╠═f25e816f-6e53-4250-a751-604f8fa463b4
# ╠═6133607f-6972-41df-8f98-9c26694f7d4a
# ╠═5d08c920-c03d-4c43-a241-428856fa7faa
