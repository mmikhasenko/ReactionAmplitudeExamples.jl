### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 2329ec14-b2eb-11ef-2df7-ad53942f462f
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add(["ThreeBodyDecays", "HadronicLineshapes",
		"Plots", "DataFrames", "Setfield"])
	# 
	using ThreeBodyDecays
	using HadronicLineshapes
	using Plots
	using LinearAlgebra
	using DataFrames
	using Setfield
end

# ╔═╡ b89d83b7-463e-416b-89f8-0b21a6c16b89
md"""
## Setup
"""

# ╔═╡ a0a1d866-59b4-4bdd-a8b8-ec58d231f26a
begin
	const mp = 0.938
	const mK = 0.548
	const mB = 5.4;
end;

# ╔═╡ 30e9aee9-5aaa-4cb6-9b94-6766c96e3900
const ms = ThreeBodyMasses(mp, mp, mK; m0=mB);

# ╔═╡ f06619b4-f60f-49bd-a833-3aec74e47036
const two_js = ThreeBodySpins(1, 1, 0; two_h0 = 0);

# ╔═╡ 5a4da2b7-ca29-469c-add6-88d4a260e74e
const tbs = ThreeBodySystem(ms, two_js);

# ╔═╡ a4a02519-f669-4541-a29a-b5ec5951de73
md"""
## Wave set
"""

# ╔═╡ 9be908cc-886c-49ee-bbfe-e29fb149ba25
df = let
	S_wave, P_wave = x2.((0, 1))	
	[
	(two_j = x2(1), two_l = S_wave, two_s = x2(1), c=cis(2π * rand())),
	(two_j = x2(0), two_l = S_wave, two_s = x2(0), c=cis(2π * rand())),
	(two_j = x2(0), two_l = P_wave, two_s = x2(1), c=cis(2π * rand())),
	(two_j = x2(1), two_l = P_wave, two_s = x2(1), c=cis(2π * rand())),
	(two_j = x2(1), two_l = P_wave, two_s = x2(0), c=cis(2π * rand())),
]
end |> DataFrame

# ╔═╡ e8226b38-b270-44d8-90ed-ff0fd57980af
# create and add decay chain object 
df_with_dc = transform(df, [:two_l, :two_s, :two_j] => ByRow() do two_l, two_s, two_j
	decay_chain = DecayChain(;
		tbs,
		k=3,
		two_j,
		Xlineshape=identity,
		Hij=RecouplingLS((two_l, two_s)), # ψ → p pbar
		HRk=RecouplingLS((two_j, two_j))  # B → ψ K
	)
	wave = ['S', 'P'][div(two_l,2)+1]
	name = "$(wave)-wave($(div(two_s,2))→$(div(two_j,2)))"
	(; name, decay_chain)
end => AsTable);

# ╔═╡ 3aa1e074-b0bc-411f-9a56-37c69c887063
model = ThreeBodyDecay(
	df_with_dc.name .=> zip(df_with_dc.c, df_with_dc.decay_chain));

# ╔═╡ 7a274027-a8f5-4aed-9d2e-6a53c808fa68
md"""
## Testing interference
"""

# ╔═╡ b0c51715-5914-4e0a-9432-87c94c6ee3ee
σs_test = randomPoint(ms);

# ╔═╡ 0f9ba6fc-6b7c-41c6-87dd-d06605874891
unpolarized_intensity(model, σs_test)

# ╔═╡ f0b86b9e-c2b1-4d90-a3d4-4338919bbf78
function interference_matrix(model, σs; objective=unpolarized_intensity)
	n = length(model)
	I_matrix = Matrix{Bool}(I, (n,n))
	intensity_m = map(Iterators.product(1:n,1:n)) do (i,j)
		_model = @set model.couplings = 
			ThreeBodyDecays.SVector(
				model.couplings .* (I_matrix[i,:] .|| I_matrix[j,:]))
		objective(_model, σs)
	end
	# leave interference: (a+b)^2-a^2-b^2
	d = diag(intensity_m)
	o = ones(Float64, n)
	x = o * d' .+ d * o' - 2diagm(d)
	D = Diagonal(intensity_m)
	D + (intensity_m - x - D) ./ 2
end

# ╔═╡ d84dad0d-36ae-4c30-a3df-56b3cb1e9016
md"""
### Unpolarized intensity (sum over proton helicity)
"""

# ╔═╡ 95de051f-45cc-4349-87e8-678e2b1a1c1d
I0 = interference_matrix(model, σs_test; objective=unpolarized_intensity);

# ╔═╡ 4d5b1fc8-3067-4f2c-9bf8-2a0ab36c6d4c
# seems vanishing for `unpolarized_intensity`
round.(I0; digits=3)

# ╔═╡ 0b791727-d382-472b-9169-4318d8f78376
begin
	_max = maximum(I0)
	clim = (-_max,_max)
	heatmap(I0; colorbar=true, clim, c=:delta, size=(600,450),
		ticks=(1:size(df_with_dc,1), df_with_dc.name))
end

# ╔═╡ 0acabb47-a24d-4e57-82ca-c4333c5b187a
md"""
### In helicity components
"""

# ╔═╡ c29ec456-43c0-4905-8666-5f83f31d8489
function amplitude_related(which_helicity_to_sum)
	f(model, σs) = sum(which_helicity_to_sum) do (l1, l2)
		abs2(amplitude(model, σs, ThreeBodySpins( l1, l2, 0; two_h0=0)))
	end
	return f
end

# ╔═╡ 4b7ffc35-eb15-4a0a-9a2e-cc4b56578862
const all_possible_helicities = Iterators.product((-1,1), (-1,1)) |> collect |> vec

# ╔═╡ 61c16841-ff86-4457-b469-5b441342722d
let # validate if matrix is computed correctly
	f = amplitude_related(all_possible_helicities)
	@assert sum(interference_matrix(model, σs_test; objective=f)) ≈ f(model, σs_test)
	@assert interference_matrix(model, σs_test; objective=f) ≈ interference_matrix(model, σs_test)
end

# ╔═╡ 7180f615-2848-4f1f-82f3-06eabad0e9f3
all_inverference_matrices = map(all_possible_helicities) do helicities
	m = interference_matrix(model, σs_test; objective=amplitude_related([helicities])) 
	round.(m; digits=3)
end

# ╔═╡ 2ee704d4-f2cf-48df-a89c-87c3ba03bb48
sum(all_inverference_matrices)

# ╔═╡ 8f1dbbbd-031b-43ce-9481-c9e870e1333e
let
	_max = maximum(abs, hcat(all_inverference_matrices...) |> vec)
	clim = (-_max, _max)
	plot(
		size=(800,160), layout=grid(1,4, widths=(0.23,0.21,0.21,0.35)),
		map(zip(all_possible_helicities, all_inverference_matrices)) do ((l1,l2),h)
			heatmap(h; colorbar=false, clim, c=:delta, title="$l1,$l2")
		end...)
	plot!(sp=1, yticks=(1:size(df_with_dc,1), df_with_dc.name))
	plot!(sp=4, colorbar=true)
end

# ╔═╡ Cell order:
# ╠═2329ec14-b2eb-11ef-2df7-ad53942f462f
# ╟─b89d83b7-463e-416b-89f8-0b21a6c16b89
# ╠═a0a1d866-59b4-4bdd-a8b8-ec58d231f26a
# ╠═30e9aee9-5aaa-4cb6-9b94-6766c96e3900
# ╠═f06619b4-f60f-49bd-a833-3aec74e47036
# ╠═5a4da2b7-ca29-469c-add6-88d4a260e74e
# ╟─a4a02519-f669-4541-a29a-b5ec5951de73
# ╠═9be908cc-886c-49ee-bbfe-e29fb149ba25
# ╠═e8226b38-b270-44d8-90ed-ff0fd57980af
# ╠═3aa1e074-b0bc-411f-9a56-37c69c887063
# ╟─7a274027-a8f5-4aed-9d2e-6a53c808fa68
# ╠═b0c51715-5914-4e0a-9432-87c94c6ee3ee
# ╠═0f9ba6fc-6b7c-41c6-87dd-d06605874891
# ╠═f0b86b9e-c2b1-4d90-a3d4-4338919bbf78
# ╟─d84dad0d-36ae-4c30-a3df-56b3cb1e9016
# ╠═95de051f-45cc-4349-87e8-678e2b1a1c1d
# ╠═4d5b1fc8-3067-4f2c-9bf8-2a0ab36c6d4c
# ╠═0b791727-d382-472b-9169-4318d8f78376
# ╟─0acabb47-a24d-4e57-82ca-c4333c5b187a
# ╠═c29ec456-43c0-4905-8666-5f83f31d8489
# ╠═4b7ffc35-eb15-4a0a-9a2e-cc4b56578862
# ╠═61c16841-ff86-4457-b469-5b441342722d
# ╠═7180f615-2848-4f1f-82f3-06eabad0e9f3
# ╠═2ee704d4-f2cf-48df-a89c-87c3ba03bb48
# ╠═8f1dbbbd-031b-43ce-9481-c9e870e1333e
