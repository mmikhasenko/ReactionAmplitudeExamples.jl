### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 7d78a382-86b0-4e45-8cf5-582620c1b131
# ╠═╡ show_logs = false
begin
	cd(mktempdir())
	using Pkg
	Pkg.activate(".")
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"),
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("PlutoUI")
	])
	# 
	using ThreeBodyDecay
	using Plots
	using PlutoUI
	import Plots.PlotMeasures: mm
end

# ╔═╡ 5b9c5979-93c4-4af4-be9f-68089147ee82
theme(:wong, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright, lab="")

# ╔═╡ 4205b30f-78c5-4e9e-9c0f-c94a10f4565f
ms = ThreeBodyMasses(0.3, 0.3, 0.3; m0=1.1)

# ╔═╡ b1e8490b-d408-4af8-a446-6cbdf1e7ba89
borderline = border(ms);

# ╔═╡ 59aa9138-16eb-4343-83e3-840bd261192c
@bind σ1 Slider(range((ms.m2+ms.m3)^2, (ms.m0-ms.m1)^2, 500),
	default=0.5, show_value=true)

# ╔═╡ ce983cb8-1d0c-423f-b040-d84bcddf658a
@bind σ2 Slider(range((ms.m3+ms.m1)^2, (ms.m0-ms.m2)^2, 500),
	default=0.5, show_value=true)

# ╔═╡ 47f81cf4-94d1-4215-81df-948fe5dfe824
let
	σs = Invariants(ms; σ1, σ2)
    (σs,α) = Kibble(σs,ms) < 0 ? (σs,1.0) : (Invariants(ms,σ1=0.7^2,σ2=0.7^2), 0.0)
	# 
	pR, p1 = four_vectors_in_binary_decay(1, 0; m1sq=σ1, m2sq=ms.m1^2, m0sq=ms.m0^2)
	p2, p3 = four_vectors_in_binary_decay(pR, cosθ23(σs, ms^2), 0;
		m1sq=ms.m2^2, m2sq=ms.m3^2)
    #
    plot(layout=grid(1,2), size=(900,350))
    plot!(sp=2, xaxis=false, yaxis=false, aspect_ratio=1)
    #
    plot!(sp=2, [0, p1[3]+1im*p1[1]], lab="", arrow=true, lw=3, α=α)
    plot!(sp=2, [0, p2[3]+1im*p2[1]], lab="", arrow=true, lw=3, α=α, seriescolor=2)
    plot!(sp=2, [0, p3[3]+1im*p3[1]], lab="", arrow=true, lw=3, α=α, seriescolor=3)
	plot!(sp=2, xlab="", ylab="")
    #
    scatter!(sp=1, [σ1], [σ2], lab="", ms=10, m=:d)
	plot!(sp=1, lab="", lc=:black, lw=3,
		getproperty.(borderline, :σ1),
		getproperty.(borderline, :σ2))
	plot!(sp=1, xlab="m23²", ylab="m12²", margins=6mm)
end

# ╔═╡ Cell order:
# ╠═7d78a382-86b0-4e45-8cf5-582620c1b131
# ╠═5b9c5979-93c4-4af4-be9f-68089147ee82
# ╠═4205b30f-78c5-4e9e-9c0f-c94a10f4565f
# ╠═b1e8490b-d408-4af8-a446-6cbdf1e7ba89
# ╠═59aa9138-16eb-4343-83e3-840bd261192c
# ╠═ce983cb8-1d0c-423f-b040-d84bcddf658a
# ╠═47f81cf4-94d1-4215-81df-948fe5dfe824
