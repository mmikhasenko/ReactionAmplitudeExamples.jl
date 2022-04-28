### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ f67af235-0ad0-421e-ae45-abf1cea4bdff
begin
	import Pkg
	Pkg.add([
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("LaTeXStrings"),
		Pkg.PackageSpec("Corpuscles"),
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl")
	])
	using Plots
	import Plots.PlotMeasures: mm
	using LaTeXStrings	
	using ThreeBodyDecay
	using Corpuscles
	using Corpuscles.Unitful
end

# ╔═╡ 0b4d7c7d-6e14-4683-95e1-456b4f4c9378
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    lw=1.2, lab="", colorbar=false)

# ╔═╡ b347c97f-facf-4281-be1b-c3ebc148bafd
begin
	const mΛc = convert(Float64,
		Particle("Lambda(c)").mass.value  / 1u"GeV*c^-2")
	const mp = convert(Float64,
		Particle("p").mass.value  / 1u"GeV*c^-2")
	const mπ = convert(Float64,
		Particle("pi").mass.value  / 1u"GeV*c^-2")
	const mK = convert(Float64,
		Particle("K").mass.value  / 1u"GeV*c^-2")
end ;

# ╔═╡ fbadb547-b5d7-4f4e-af0a-2af068c15e6a
const ms = ThreeBodyMasses(mp, mK, mπ; m0=mΛc)

# ╔═╡ c2d8a53b-50f9-40af-a33c-96b1aeef24cf
function field((σ1,σ3); wr=wr(3,1,0))
	σs = Invariants(ms; σ3,σ1)
	Kibble(σs, ms^2) > 0 && return (NaN,NaN)
	c = cosζ(wr, σs, ms^2)
	(c, sqrt(1-c^2) * (ispositive(wr) ? 1 : -1))
end

# ╔═╡ 4b4d6a4b-862a-40cb-980b-6e604638b29e
begin
	function unscale((x,y))
		lx = lims1(ms)
		ly = lims3(ms)
		(lx[1] + (lx[2]-lx[1])*x,
		 ly[1] + (ly[2]-ly[1])*y)
	end
	function scale((u,v))
		lx = lims1(ms)
		ly = lims3(ms)
		((u-lx[1])/(lx[2]-lx[1]),
		 (v-ly[1])/(ly[2]-ly[1]))
	end
end

# ╔═╡ 4d1a05dd-3477-401f-94f4-e2f8077d0603
begin
	function transform2tuple(x,yv)
		y1,y2 = yv[:,1],yv[:,2]
		vcat(Tuple.(zip(x,y1)), Tuple.(zip(reverse(x),reverse(y2))))
	end
	swap((x,y)) = (y,x)
end

# ╔═╡ 08af1f39-dcfe-4a19-b408-827c32ba7c4f
let length=18
	xv = range(0,1; length)
	yv = range(0,1; length)
	xyv = Iterators.product(xv, yv)
	plot(layout=grid(1,2), size=(700,350),
		title=[L"\alpha(\Delta)" L"\alpha(\Lambda)"],
		xlab=L"\sim m(\pi K)", ylab=L"\sim m(pK)",
		xaxis=nothing, yaxis=nothing, bottom_margin=3mm)
	# 
	quiver!(sp=1, getindex.(xyv,1), getindex.(xyv,2),
		quiver=(x,y)->field(unscale((x,y)); wr=wr(3,1,0)) ./ length,
		c=2, arrow=(:closed,), aspectratio=1)
	plot!(sp=1, map(scale, transform2tuple(border31(ms)...)), lc=:black)
	# 
	quiver!(sp=2, getindex.(xyv,1), getindex.(xyv,2),
		quiver=(x,y)->field(unscale((x,y)); wr=wr(2,1,0)) ./ length,
		c=3, arrow=(:closed,), aspectratio=1)
	plot!(sp=2, map(scale, transform2tuple(border31(ms)...)), lc=:black)
	# 
end

# ╔═╡ Cell order:
# ╠═f67af235-0ad0-421e-ae45-abf1cea4bdff
# ╠═0b4d7c7d-6e14-4683-95e1-456b4f4c9378
# ╠═b347c97f-facf-4281-be1b-c3ebc148bafd
# ╠═fbadb547-b5d7-4f4e-af0a-2af068c15e6a
# ╠═c2d8a53b-50f9-40af-a33c-96b1aeef24cf
# ╟─4b4d6a4b-862a-40cb-980b-6e604638b29e
# ╟─4d1a05dd-3477-401f-94f4-e2f8077d0603
# ╠═08af1f39-dcfe-4a19-b408-827c32ba7c4f
