### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 90c57cae-40c5-11ed-3374-45e9dd35d5de
# ╠═╡ show_logs = false
begin
	import Pkg
	Pkg.add(
		[Pkg.PackageSpec(
			url="https://github.com/mmikhasenko/EffectiveRangeExpansion.jl"),
		Pkg.PackageSpec("QuadGK"),
		Pkg.PackageSpec("Plots")])
	
	using EffectiveRangeExpansion
	using QuadGK
	using Plots
end

# ╔═╡ a3bf6dcb-d6c3-44f1-b90b-b23f66c5a75a
md"""
# Effective range of the two-body loop 
"""

# ╔═╡ d2c19af4-724a-48ab-9a33-305ad36fd7e5
arg(z) = atan(imag(z), real(z))

# ╔═╡ 82773afd-8fde-41fb-b50b-d811fbf10eab
begin
	const m₁ = 1
	const m₂ = 1
end

# ╔═╡ ff6cf97f-674b-431b-b3fa-b02f1841194c
begin
	const eth = m₁+m₂
	const sth = eth^2
end

# ╔═╡ e4ced65c-9eb6-4c15-aaf1-dcc9cf6bf9c5
md"""
The direction of the branch cut in the phase space function is adjusted by the phase factor $e^{-i\phi}$. The phase $\phi=0$ corresponds to the cut to the right,
$\phi=-\pi/2$ puts the cut straight down
"""

# ╔═╡ 31797bb4-2bb0-4525-8993-ac51d6641ceb
const ϕ = 0

# ╔═╡ 8a7077b4-eb27-43e7-9c0f-98c9aa0be48e
ρ(e) = 1im*cis(ϕ/2) * sqrt(cis(-ϕ)*((m₁+m₂)-e)) *
	sqrt((m₁+m₂)+e) * 
	sqrt(e-(m₁-m₂)) *
	sqrt(e+(m₁-m₂)) / e^2

# ╔═╡ d57d4e8b-3e76-4d60-88d5-99c215697dca
md"""
To make sure that the $\rho$ is defined with correct sign, it is checked against the simpler equation.
"""

# ╔═╡ 7e925c5b-520b-4dba-982f-1aaa16490ad7
ρR(e) = 1im*sqrt((m₁+m₂)^2-e^2)/e

# ╔═╡ 1d09e3eb-d541-40af-b3e7-a38bbdf3d67b
begin
	@assert ρR(1+1im) ≈ ρ(1+1im)
	@assert ρR(3+1im) ≈ ρ(3+1im)
end

# ╔═╡ 934ce5cc-6033-48a6-89ce-5b576acc9de1
heatmap(-0.1:0.01:3, -1:0.01:1, (x,y)->arg(1im*ρ(x+1im*y)), c=:twilight)

# ╔═╡ 72a252a4-ba7d-42ad-b730-262aab8d0cd7
md"""
The Chew-Mandestam function is defined by the dispersion integral,

$C(s) = \frac{s}{\pi}\int_{s_\text{th}}^\infty\frac{\rho(s')}{s'-s-i0} \mathrm{d}s'\,,$

with $\text{Im}\,C(s) = \rho(s)$
"""

# ╔═╡ 760cc78c-1d4f-4a7e-9256-84c4ce23a176
C₀(e) = e^2/π * quadgk(s′->ρ(sqrt(s′))/s′/(s′-e^2-1e-7im), sth, Inf)[1]

# ╔═╡ 3fc670db-7777-4a04-8c29-74870b872b53
md"""
The subtraction point is shifted to the threshold by redefining the function.
"""

# ╔═╡ 63ad7fd8-0644-4a4d-a12a-e3b9d7243836
const c₀ = real(C₀(eth))

# ╔═╡ 3afbbcd6-4be4-486a-b09c-298b4fbbdffc
C(s) = C₀(s)-c₀

# ╔═╡ e269df3b-d8c7-4a0c-916c-1a756562ed52
md"""
The analytic struction is determined by the integral form.
It is a cut starting from the threshold stright to the right along the real axis.
"""

# ╔═╡ 7b820569-580e-4ae5-8c58-60a72fdd1711
heatmap(-0.1:0.01:3, -1:0.01:1, (x,y)->arg(C(x+1im*y)), c=:twilight)

# ╔═╡ 74be6176-f715-414f-a1d8-d5a954ae4825
begin
	plot(grid=false, frame=:box, xlabel="√s", ylab="Amplitude")
	plot!(e->imag(C(e)), 1, 3, l=:red, lab="Im")
	plot!(e->real(C(e)), 1, 3, l=:black, lab="Re")
end

# ╔═╡ a143f531-4c00-4464-89cf-619fc1ad4e14
md"""
The scattering parameters are defined by the expansion,

$f(s) = N (a⁻¹ + \frac{r k(s)^2}{2}  - i k(s))$

The values are computed numerically using the technique of the Cauchy integrals, see `[EffectiveRangeExpansion.jl]`
"""

# ╔═╡ f1265120-b223-4796-8a6e-0d4f71b1491a
md"""
First, the test shows that for $f = -ik$, the scattering parameters are zero
"""

# ╔═╡ 6950bb68-1990-484c-9e05-168f697f2875
ere₀ = let
	f(x) = -1im*ρ(x+eth)
	k(x) = ρR(x+eth)
	effectiverangeexpansion(f, k,
		ComplexBranchPointExpansion(CircularIntegral(0.01)))
end

# ╔═╡ c85daf45-9cec-40c2-985a-da3d7aad8d32
md"""
The effective range is **positive** the function $f = -C(s)$.
It is the standard denominator of the Breit-Wigner with running mass,

$m^2-s- im\Gamma(s)\,\to\,m^2-s- ig^2 \rho(s)\,\to\,m^2-s - g^2C(s)$,

where the potential terms, $g^2/(m^2-s)$ is omitted.
"""

# ╔═╡ e732e2ba-3c26-4be5-97cf-a33b11818fff
ere = let
	f(x) = -C(x+eth)
	k(x) = ρR(x+eth)
	effectiverangeexpansion(f, k,
		ComplexBranchPointExpansion(CircularIntegral(0.01)))
end

# ╔═╡ Cell order:
# ╟─a3bf6dcb-d6c3-44f1-b90b-b23f66c5a75a
# ╠═90c57cae-40c5-11ed-3374-45e9dd35d5de
# ╠═d2c19af4-724a-48ab-9a33-305ad36fd7e5
# ╠═82773afd-8fde-41fb-b50b-d811fbf10eab
# ╠═ff6cf97f-674b-431b-b3fa-b02f1841194c
# ╟─e4ced65c-9eb6-4c15-aaf1-dcc9cf6bf9c5
# ╠═31797bb4-2bb0-4525-8993-ac51d6641ceb
# ╠═8a7077b4-eb27-43e7-9c0f-98c9aa0be48e
# ╟─d57d4e8b-3e76-4d60-88d5-99c215697dca
# ╠═7e925c5b-520b-4dba-982f-1aaa16490ad7
# ╠═1d09e3eb-d541-40af-b3e7-a38bbdf3d67b
# ╠═934ce5cc-6033-48a6-89ce-5b576acc9de1
# ╟─72a252a4-ba7d-42ad-b730-262aab8d0cd7
# ╠═760cc78c-1d4f-4a7e-9256-84c4ce23a176
# ╟─3fc670db-7777-4a04-8c29-74870b872b53
# ╠═63ad7fd8-0644-4a4d-a12a-e3b9d7243836
# ╠═3afbbcd6-4be4-486a-b09c-298b4fbbdffc
# ╟─e269df3b-d8c7-4a0c-916c-1a756562ed52
# ╠═7b820569-580e-4ae5-8c58-60a72fdd1711
# ╠═74be6176-f715-414f-a1d8-d5a954ae4825
# ╟─a143f531-4c00-4464-89cf-619fc1ad4e14
# ╟─f1265120-b223-4796-8a6e-0d4f71b1491a
# ╠═6950bb68-1990-484c-9e05-168f697f2875
# ╟─c85daf45-9cec-40c2-985a-da3d7aad8d32
# ╠═e732e2ba-3c26-4be5-97cf-a33b11818fff
