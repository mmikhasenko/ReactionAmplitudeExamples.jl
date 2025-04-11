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

# ╔═╡ b077d8ea-1620-11f0-2d73-ad9dda501211
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ScatteringKMatrix.jl"),
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("PlutoUI"),
		Pkg.PackageSpec("Parameters")
	])
	# 
	using ScatteringKMatrix
	using ScatteringKMatrix.StaticArrays
	using Plots
	using Parameters
	using PlutoUI
end

# ╔═╡ 710bd318-025c-4761-bc8b-51c07d1dd40a
theme(:boxed)

# ╔═╡ f713eddb-035b-4f97-8f37-db147ad111fd
channels = SVector(
        TwoBodyChannel(1.1, 1.1),
        TwoBodyChannel(2.2, 2.2),
    )

# ╔═╡ 2d90bc7e-0e2d-4496-8209-90a2c7637299
function build_model(pars)
	@unpack M1,g1,g2 = pars
	@unpack M2,h1,h2 = pars
	@unpack α1, α2 = pars
	MG = [
		(M = M1, gs = [g1, g2]),
		(M = M2, gs = [h1, h2])
	]
	K = KMatrix(MG)
	T = TMatrix(K, channels)
	# 
	ProductionAmplitude(T, [α1, α2])
end

# ╔═╡ d4927b6a-3fb3-4709-a785-f41ad1faae5f
md"""

### First K-matrix pole
- M1 = $(@bind M1 Slider(5.0:0.01:6.0, show_value=true, default=5.3))
- g1 = $(@bind g1 Slider(1.0:0.01:3.0, show_value=true, default=2.2))
- g2 = $(@bind g2 Slider(1.0:0.01:3.0, show_value=true, default=1.48))

### Second K-matrix pole
- M2 = $(@bind M2 Slider(7.0:0.01:9.0, show_value=true, default=8.3))
- h1 = $(@bind h1 Slider(2.0:0.01:4.0, show_value=true, default=3.2))
- h2 = $(@bind h2 Slider(2.0:0.01:4.0, show_value=true, default=3.48))

### Production
- α1 = $(@bind α1 Slider(0.5:0.01:2.0, show_value=true, default=1.1))
- α2\_abs = $(@bind α2_abs Slider(0.1:0.01:2.0, show_value=true, default=1.3))
- α2\_phi = $(@bind α2_phi Slider(-π:0.01:-π, show_value=true, default=1.3))

"""

# ╔═╡ 6ef38e97-9b55-4c1d-8a0d-bfbc937370c4
F = let
	pars = (;
		M1, g1, g2,
		M2, h1, h2,
		α1,
	 	α2=α2_abs*cis(α2_phi))
	build_model(pars)
end

# ╔═╡ 75e9a59c-a537-477e-aef9-6581acdf1e8d
A(m) = amplitude(F, m)

# ╔═╡ 5eca88a9-f724-4420-a38a-bf849df20ae9
mth, m_max = 4.4, 10

# ╔═╡ 6104256c-7df7-42c5-ac3a-386e0f337c4e
begin
	plot()
	plot!(m->abs2(A(m)[1]), mth, m_max)
	plot!(m->abs2(production_pole(F, m, 1)[1]), mth, m_max, fill=0, fillalpha=0.3)
	plot!(m->abs2(production_pole(F, m, 2)[1]), mth, m_max, fill=0, fillalpha=0.3)
end

# ╔═╡ Cell order:
# ╠═b077d8ea-1620-11f0-2d73-ad9dda501211
# ╠═710bd318-025c-4761-bc8b-51c07d1dd40a
# ╠═f713eddb-035b-4f97-8f37-db147ad111fd
# ╠═2d90bc7e-0e2d-4496-8209-90a2c7637299
# ╟─d4927b6a-3fb3-4709-a785-f41ad1faae5f
# ╟─6ef38e97-9b55-4c1d-8a0d-bfbc937370c4
# ╠═75e9a59c-a537-477e-aef9-6581acdf1e8d
# ╠═5eca88a9-f724-4420-a38a-bf849df20ae9
# ╠═6104256c-7df7-42c5-ac3a-386e0f337c4e
