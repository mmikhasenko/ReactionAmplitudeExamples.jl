### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ ce0385b0-b12f-11ed-080d-1fe672667cfc
begin
	using Pkg
	cd("C://Users//mikha//Documents//X2DDpi")
	Pkg.activate(".")
end

# ╔═╡ c01e2ec1-4785-45cc-9972-457b261d5316
using Plots

# ╔═╡ bd0fc88d-655f-4545-83c7-baede8be1fe9
using AlgebraPDF

# ╔═╡ 60d8993c-31ca-4483-bd6f-c46fa3f8323a
theme(:wong2, frame=:box, grid=:false,
	xlims=(:auto,:auto), ylims=(0,:auto), lw=2)

# ╔═╡ 7ca7226f-60df-4180-ab10-fdf1f4d660c7
md"""
## Color wave function
"""

# ╔═╡ 88766e27-a6ec-4ebd-ac38-4ee2b6b3f1dc
ρ_bb_col = let
	ρ_bb_3 = Normalized(
		FunctionWithParameters((x;p)->x^2*exp(-x/0.08); p=∅),
		(0,2))
	ρ_bb_6 = Normalized(
		FunctionWithParameters((x;p)->x^2*exp(-x/0.2); p=∅),
		(0,2))
	(c3=1,)*ρ_bb_3 + (c6=1/30,)*ρ_bb_6
end

# ╔═╡ f5071865-ddeb-43a8-9da7-02aca0f41bff
let
	plot(ρ_bb_col, lab="total")
	plot!(ρ_bb_col[1], lab="triplet")
	plot!(ρ_bb_col[2], norm=30, lab="sixtet")
end

# ╔═╡ 23ebf3e8-a0ea-40ad-b088-9cd15d6ed5c2
begin
	N = 5000
	N3 = round(Int, N * ρ_bb_col.αs[1])
	N6 = round(Int, N * ρ_bb_col.αs[2])
	# 
	d3 = rand(ρ_bb_col.fs[1], N3)
	d6 = rand(ρ_bb_col.fs[2], N6)
end

# ╔═╡ 9db090e7-e42c-4052-8ced-bc129783e741
begin
	plot(xlims=(-1,1), ylims=(-1,1), size=(300,300))
	scatter!(d3 .* cis.( 2π .* rand(length(d3))), frame=:origin, c=2, α=0.3, lab="triplet (3)", ms=3)
	scatter!(d6 .* cis.( 2π .* rand(length(d6))), frame=:origin, c=3, α=0.7, lab="sixtet (6)", ms=3)
	plot!(xlab="x", ylab="y", aspect_ratio=1)
	savefig("TQQ_color_scatter.pdf")
	plot!()
end

# ╔═╡ a9098d3d-58df-4603-add2-93e29bb83f65
md"""
## Density wave function
"""

# ╔═╡ 41d3ac4a-e3e5-4f25-893a-a351765e506b
begin
	ρ_bb_spa = Normalized(
		FunctionWithParameters((x;p)->x^2*exp(-x/0.08); p=∅),
		(0,2))
	ρ_bq_spa = Normalized(
		FunctionWithParameters((x;p)->x^2*exp(-x^2/0.4); p=∅),
		(0,2))
end

# ╔═╡ 7e6a40c1-f4d0-4763-b257-f1b59fbf3cf9
let
	plot()
	plot!(ρ_bb_spa, lab="bb")
	plot!(ρ_bq_spa, norm=30, lab="bq")
end

# ╔═╡ 054eabe1-31c9-4d12-9c63-4061d5a2afca
begin
	N_spa = 1000
	# 
	d_bb = rand(ρ_bb_spa, N_spa)
	d_bq = rand(ρ_bq_spa, N_spa)
end;

# ╔═╡ 69d3ffa4-bdd6-40b0-86d1-731e8a0e41e7
begin
	plot(xlims=(-1.5, 1.5), ylims=(-1.5, 1.5), size=(300,300))
	scatter!(d_bb .* cis.( 2π .* rand(length(d_bb))), frame=:origin, c=2, α=0.5, lab="bb", ms=3)
	scatter!(d_bq .* cis.( 2π .* rand(length(d_bq))), frame=:origin, c=3, α=0.5, lab="qq", ms=3)
	plot!(xlab="x", ylab="y", aspect_ratio=1)
	savefig("TQQ_spa_scatter.pdf")
	plot!()
end

# ╔═╡ Cell order:
# ╠═ce0385b0-b12f-11ed-080d-1fe672667cfc
# ╠═60d8993c-31ca-4483-bd6f-c46fa3f8323a
# ╠═c01e2ec1-4785-45cc-9972-457b261d5316
# ╠═bd0fc88d-655f-4545-83c7-baede8be1fe9
# ╟─7ca7226f-60df-4180-ab10-fdf1f4d660c7
# ╠═88766e27-a6ec-4ebd-ac38-4ee2b6b3f1dc
# ╠═f5071865-ddeb-43a8-9da7-02aca0f41bff
# ╠═23ebf3e8-a0ea-40ad-b088-9cd15d6ed5c2
# ╠═9db090e7-e42c-4052-8ced-bc129783e741
# ╟─a9098d3d-58df-4603-add2-93e29bb83f65
# ╠═41d3ac4a-e3e5-4f25-893a-a351765e506b
# ╠═7e6a40c1-f4d0-4763-b257-f1b59fbf3cf9
# ╠═054eabe1-31c9-4d12-9c63-4061d5a2afca
# ╠═69d3ffa4-bdd6-40b0-86d1-731e8a0e41e7
