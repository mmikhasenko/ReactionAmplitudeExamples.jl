### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ e8300b76-37a6-11ef-3955-f5c442c10a6a
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(".")
	
	using ThreeBodyDecays
	using Plots
end

# ╔═╡ 01ac5398-e641-45b5-b527-86aaf8553d98
theme(:wong2, lab="", grid=false, frame=:box)

# ╔═╡ a6442249-dacf-4d6d-bbbc-bebcde52b5f9
const ms = ThreeBodyMasses(0.14,0.14,0.14; m0=1.4)

# ╔═╡ 20fa77d7-8b3e-4a1d-9901-96609abcdd50
begin
	data = y2σs.(eachslice(rand(3,1_000_000), dims=2), Ref(ms))
	filter!(Base.Fix2(isphysical, ms), data)
end;

# ╔═╡ 20de7e99-a96a-437b-80e3-b2d27032c010
begin
	const mω = 0.789
	const Δω = 0.03
end

# ╔═╡ 3dec26f9-7ad9-45c8-a3c0-81effde0f2be
reduced_data = filter(data) do (σ1,σ2,σ3)
	f1 = abs(sqrt(σ1)-mω) < Δω
	f2 = abs(sqrt(σ2)-mω) < Δω
	(f1 || f2) && (f1 != f2)
end;

# ╔═╡ c83ec623-afcb-4d65-a098-bfcc9263f18d
begin
	histogram2d(reduced_data, bins=150)
	plot!(border12(ms), l=(5,:black))
end

# ╔═╡ cdfaf595-f6ac-47af-8956-f63ee88e3b87
length(reduced_data) / length(data)

# ╔═╡ Cell order:
# ╠═e8300b76-37a6-11ef-3955-f5c442c10a6a
# ╠═01ac5398-e641-45b5-b527-86aaf8553d98
# ╠═a6442249-dacf-4d6d-bbbc-bebcde52b5f9
# ╠═20fa77d7-8b3e-4a1d-9901-96609abcdd50
# ╠═20de7e99-a96a-437b-80e3-b2d27032c010
# ╠═3dec26f9-7ad9-45c8-a3c0-81effde0f2be
# ╠═c83ec623-afcb-4d65-a098-bfcc9263f18d
# ╠═cdfaf595-f6ac-47af-8956-f63ee88e3b87
