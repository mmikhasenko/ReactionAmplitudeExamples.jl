### A Pluto.jl notebook ###
# v0.16.1

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

# ╔═╡ 43870f95-68c5-459b-b18c-e8571b87ed2d
begin
    import Pkg
    Pkg.add(Pkg.PackageSpec(url="https://github.com/mmikhasenko/Triangles"))
    using Plots, LaTeXStrings, TriangleSingularity, PlutoUI
end

# ╔═╡ 8d16c0a0-685f-44cf-afed-2fa3e23e4d05
md"""
# $X \to D^{*0}D^0$ triangle
"""

# ╔═╡ 33f0765d-5cd0-4c01-a0b6-c587506102e8
# 2021-11-16 Misha Mikhasenko

# ╔═╡ 6d911b63-eb21-49b1-964a-aeab2433cf9b
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto), lab="",
    xlab=L"m(D^0D^{*0})\,\,[\mathrm{MeV}]")

# ╔═╡ 64603ace-d54b-4040-95b2-0f54f1b69bc1
begin
	const mπ⁰ = 0.1349768
	const mπ⁺ = 0.13957039
	const mD⁰ = 1.86483
	const mD⁺ = 1.86965
	# 
	const mDˣ⁺ = 1.86483+145.4258e-3 # m(D) + Δm(D*,D) from PDG
	const mDˣ⁰ = 2.00685
	# 
	const ΓDˣ⁺ = 83.4e-6
	const ΓDˣ⁰ = 55.2e-6
	const mγ = 0.0
end

# ╔═╡ d20805dc-337f-4055-b2f2-4081dead3c27
tr(δe,mpi,Γ=ΓDˣ⁰) = 
	triangleloop(
		mDˣ⁰^2-1im*mDˣ⁰*Γ, mD⁰^2, mpi^2, mDˣ⁰^2, mD⁰^2, (mDˣ⁰+mD⁰+δe*1e-3)^2)

# ╔═╡ 11c03102-842f-42f1-adeb-327c2f6156eb
@bind _mπ Slider(range(0.100, 0.250,length=76); show_value=true, default=0.140)

# ╔═╡ 660cc2c4-174e-4634-9925-c521231cbfc7
begin
	plot()
	plot!(δe->imag(tr(δe,_mπ)), -1.5,1.5,
		lab=L"\Gamma_{D^{*0}}=80\,\mathrm{keV}")
	plot!(δe->imag(tr(δe,_mπ,1e-7)), -1.5,1.5, 
		lab=L"\Gamma_{D^{*0}}=0\,\mathrm{keV}")
end

# ╔═╡ Cell order:
# ╠═43870f95-68c5-459b-b18c-e8571b87ed2d
# ╟─8d16c0a0-685f-44cf-afed-2fa3e23e4d05
# ╠═33f0765d-5cd0-4c01-a0b6-c587506102e8
# ╠═6d911b63-eb21-49b1-964a-aeab2433cf9b
# ╠═64603ace-d54b-4040-95b2-0f54f1b69bc1
# ╠═d20805dc-337f-4055-b2f2-4081dead3c27
# ╠═11c03102-842f-42f1-adeb-327c2f6156eb
# ╠═660cc2c4-174e-4634-9925-c521231cbfc7
