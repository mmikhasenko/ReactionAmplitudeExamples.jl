### A Pluto.jl notebook ###
# v0.11.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ cff8dfe0-e3d3-11ea-20ce-6bf12d4cdbf0
using Plots

# ╔═╡ 454149e0-e3d4-11ea-3e4f-f10ff6a7af9c
using PlutoUI

# ╔═╡ cdf07ad0-e3d5-11ea-39bf-cbd4c210bde9
sth = 0.9;

# ╔═╡ 6d9d0e50-e3d5-11ea-3f00-d7e117d169d7
md"""Γ: $(@bind Γ Slider(range(0.1, 0.3, length=100);
	default=0.2, show_value=true))"""

# ╔═╡ be4a5ace-e3d3-11ea-004f-2fe55918802c
md"""m: $(@bind m Slider(range(0.95, 1.0, length=150);
	default=1.0, show_value=true))"""

# ╔═╡ 7dffaf40-e3d1-11ea-04fa-6b336eadb946
BW(s) = m*Γ/(m^2-s-m*Γ*sqrt(-1+sth/s))

# ╔═╡ d2f5f340-e3d3-11ea-2e43-059ba1ce3b75
let
	xv = range(0.9, 1.0, length=100)
	yv = range(-0.05, 0.05, length=100)
	calv = BW.((xv'.+1im.*yv).^2)
	contour(xv,yv,log.(abs2.(calv)), c=:viridis, colorbar=false)
	plot!([sqrt(sth), 1.0+0im], l=(2,:red), lab="cut")
	scatter!([sqrt(sth)+0im], m=(4,:red), lab="")
end

# ╔═╡ Cell order:
# ╠═cff8dfe0-e3d3-11ea-20ce-6bf12d4cdbf0
# ╠═454149e0-e3d4-11ea-3e4f-f10ff6a7af9c
# ╠═cdf07ad0-e3d5-11ea-39bf-cbd4c210bde9
# ╟─6d9d0e50-e3d5-11ea-3f00-d7e117d169d7
# ╟─be4a5ace-e3d3-11ea-004f-2fe55918802c
# ╠═7dffaf40-e3d1-11ea-04fa-6b336eadb946
# ╟─d2f5f340-e3d3-11ea-2e43-059ba1ce3b75
