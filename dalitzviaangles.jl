### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ f7e83510-9fb5-11ec-3f86-ad189b9ee734
begin
	import Pkg
	Pkg.add(Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"))
	Pkg.add("Plots")
	# 
	using ThreeBodyDecay
	using Plots
end

# ╔═╡ 7b93b085-9d13-4fd2-abb0-c39bbcb823e1
ms = ThreeBodyMasses(1.1,0.8,0.3; m0=4.2)

# ╔═╡ 63a940ea-dcd3-42dd-bdd2-847a1be17ffd
ms² = ms^2

# ╔═╡ 6c721d01-a80a-486e-a6e1-29fc9d78db84
data = flatDalitzPlotSample(ms, Nev=1000)

# ╔═╡ d97defca-748c-45b4-82ff-48c2a1eb9867
zst = [(cosθ12=cosθ12(x, ms²), cosθ23=cosθ23(x, ms²)) for x in data]

# ╔═╡ e0abcd0d-eb69-4cd9-8409-9d7c139f2ff3
scatter(getindex.(data,1), getindex.(data,2))

# ╔═╡ 0f5bb842-f86a-446e-a65a-50b45438233c
scatter(getindex.(zst,1), getindex.(zst,2))

# ╔═╡ Cell order:
# ╠═f7e83510-9fb5-11ec-3f86-ad189b9ee734
# ╠═7b93b085-9d13-4fd2-abb0-c39bbcb823e1
# ╠═63a940ea-dcd3-42dd-bdd2-847a1be17ffd
# ╠═6c721d01-a80a-486e-a6e1-29fc9d78db84
# ╠═d97defca-748c-45b4-82ff-48c2a1eb9867
# ╠═e0abcd0d-eb69-4cd9-8409-9d7c139f2ff3
# ╠═0f5bb842-f86a-446e-a65a-50b45438233c
