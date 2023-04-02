### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 1c6e06ae-5747-11ec-0f31-9347d77fe04b
begin
	import Pkg
	Pkg.add("AlgebraPDF")
	Pkg.add("QuadGK")
	Pkg.add("Plots")
	Pkg.add("LaTeXStrings")
	# 
	using AlgebraPDF
	using Plots
	using QuadGK
	using LaTeXStrings
end

# ╔═╡ d9ca632b-2ab0-4f2d-bb56-17aa616665b5
begin
	import Plots.PlotMeasures: mm
	theme(:wong2, frame=:box, grid=false, minorticks=true, 
	    guidefontvalign=:top, guidefonthalign=:right,
	    xlim=(:auto,:auto), ylim=(:auto,:auto),
	    xlab=L"m_{3\pi}", lw=1.2)
end

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
	UnitarizedBackground(B,T,ρ) = UnitarizedBackground(B,T,ρ,1e-6)
	#
	AlgebraPDF.pars(ub::UnitarizedBackground, isfree::Bool) =
		pars(ub.B, isfree) + pars(ub.T, isfree) + pars(ub.ρ, isfree)
	# 
	AlgebraPDF.func(ub::UnitarizedBackground, x::NumberOrTuple; p=pars(ub)) =
		#
		func(ub.B,x;p) + func(ub.T,x;p)*
			quadgk(x′->func(ub.B,x′;p)*ub.ρ(x′)/(x′-x-1im*ub.ϵ), 0.6^2, Inf)[1]/π
end

# ╔═╡ b1577a00-af60-47c5-936e-1bc1134a133a
begin
	@makefuntype Background(s;p) = s < sth ? 0.0 : (s-sth)^p.α*exp(-p.β*s)
	b = Background((α=5.1, β=1.2))
end

# ╔═╡ e78d6e79-2d1f-4471-8489-1d7d8ed33dd8
plot(e->b(e^2),0.6,2.5)

# ╔═╡ fbaa22eb-041a-4525-bd11-0b6ca0f0dc4f
begin
	@makefuntype Signal(s;p) = p.m₀*p.Γ₀/(p.m₀^2-s-1im*p.m₀*p.Γ₀)
	sig = Signal((m₀=1.8, Γ₀=0.15))
end

# ╔═╡ 2041e2f1-87b7-4cbe-accd-c462a71212a9
plot(e->abs2(sig(e^2)),0.6,2.5)

# ╔═╡ a4d25e58-2815-4610-b90f-2225023a4609
ρ(s) = 1.0

# ╔═╡ bb05f42f-d70b-4824-978a-8e398d164d3d
ub = UnitarizedBackground(b,sig,ρ)

# ╔═╡ 00374615-e9e5-4774-a8d9-2377416ce042
let
	plot()
	plot!(e->abs2(ub(e^2)),0.6,2.5, lab="Unitarized(B by S)", lw=2)
	plot!(e->abs2(b(e^2)),0.6,2.5, lab="Background")
	plot!(e->abs2(sig(e^2))*70,0.6,2.5, lab="Signal")
end

# ╔═╡ d5182be7-b732-4b19-95d0-de1e4bc8160c
sub = sig*(c1=70.1,) + ub*(c2=1.0,)

# ╔═╡ f25e816f-6e53-4250-a751-604f8fa463b4
c1=10

# ╔═╡ 5d08c920-c03d-4c43-a241-428856fa7faa
let
	d = updatepar(sub,:c1, c1)
	plot(e->abs2(d(e^2)),0.6,2.5, lab="S+Unitarized(B)", lw=2)
	plot!(e->abs2((sig*(;c1))(e^2)),0.6,2.5, lab="S")
end

# ╔═╡ Cell order:
# ╠═1c6e06ae-5747-11ec-0f31-9347d77fe04b
# ╠═d9ca632b-2ab0-4f2d-bb56-17aa616665b5
# ╠═18663f9a-5d04-4530-84e7-cadf474643da
# ╠═b1577a00-af60-47c5-936e-1bc1134a133a
# ╠═e78d6e79-2d1f-4471-8489-1d7d8ed33dd8
# ╠═fbaa22eb-041a-4525-bd11-0b6ca0f0dc4f
# ╠═2041e2f1-87b7-4cbe-accd-c462a71212a9
# ╟─ab563532-193c-4cdf-a837-13124ce4da78
# ╠═2fdf557e-d2bb-4739-9004-ced593298f09
# ╠═a4d25e58-2815-4610-b90f-2225023a4609
# ╟─bb05f42f-d70b-4824-978a-8e398d164d3d
# ╠═00374615-e9e5-4774-a8d9-2377416ce042
# ╠═d5182be7-b732-4b19-95d0-de1e4bc8160c
# ╠═f25e816f-6e53-4250-a751-604f8fa463b4
# ╠═5d08c920-c03d-4c43-a241-428856fa7faa
