### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ 3d7df13c-75a7-46ac-b28e-18f263de6af7
begin
	using Pkg
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"),
		Pkg.PackageSpec("PartialWaveFunctions"),
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("PlutoUI")
		])
	
	using ThreeBodyDecay
	using PartialWaveFunctions
	using Plots
	using PlutoUI
end

# ╔═╡ 8dba788b-f74f-441c-a712-6116761a4a83
md"""
### Angular distribution for $J/\psi\to (\pi\pi)_L\,\pi\,L$-wave
"""

# ╔═╡ 66d3aa29-51c8-49bd-8df3-bc6469d18901
theme(:wong2, minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))

# ╔═╡ 54f6f89e-f496-4ad9-8047-44e9dbd43ef4
A(z,s,λ) = (2s+1)*wignerd(s,λ,0,z)*clebschgordan(s,0,s,λ,1,λ)

# ╔═╡ a8c0f283-92fe-4433-b5dc-341b14c7a0a2
Z(z,s) = sum(abs2, A(z,s,λ) for λ in -1:1)

# ╔═╡ 9dcf87e4-beff-4b66-9118-6b6656e35e46
begin
	plot(xlab="cosθ", title="angular distribution J/ψ → (ππ)ₛ π s-wave")
	plot!(z->Z(z,1), -1, 1, lab="P-wave")
	plot!(z->Z(z,3), -1, 1, lab="F-wave")
end

# ╔═╡ 876e7087-3948-4f84-9dcd-a4156c563983
md"""
#### Dalitz Plot for different waves
I would like to see how a general patter differs between P-wave and F-wave.
For the quick check:
 * I neglect the interference,
 * Consider some broad resonance $m=1.5\,$GeV, $m=0.3\,$GeV.
"""

# ╔═╡ 86ba31c9-57bf-4cd2-8eb0-d66be0ddc89f
ms = ThreeBodyMasses(0.15,0.15,0.15; m0=3.09)

# ╔═╡ e3206b0d-277a-44e8-a26b-355aeb1dcae2
function I(σs; p)
	abs2(BW(σs.σ1, p.m, p.Γ))*Z(cosθ23(σs, ms^2), p.s)+
		abs2(BW(σs.σ2, p.m, p.Γ))*Z(cosθ31(σs, ms^2), p.s)+
		abs2(BW(σs.σ3, p.m, p.Γ))*Z(cosθ12(σs, ms^2), p.s)
end

# ╔═╡ b60ec657-a8bd-4229-9848-4f3a7fab7fde
@bind mR Slider(0:0.1:3, default=1.5, show_value=true)

# ╔═╡ fc53f25c-c6d7-402c-9d1f-2f41daed0fe1
@bind ΓR Slider(0.05:0.01:0.5, default=0.2, show_value=true)

# ╔═╡ 69a05a98-f0e6-49de-be57-42b51dd0b551
begin
	plot(layout=grid(1,2), size=(650,300), grid=false, frame=:box)
	# 	
	plot!(sp=1, σs->I(σs; p=(s=1, m=mR, Γ=ΓR)), ms, c=:viridis)
	plot!(sp=1, border12(ms), l=(2,:black), lab="", title="P-wave")
	# 	
	plot!(sp=2, σs->I(σs; p=(s=3, m=mR, Γ=ΓR)), ms, c=:viridis)
	plot!(sp=2, border12(ms), l=(2,:black), lab="", title="F-wave")
end

# ╔═╡ baa08139-93f6-4c6f-9fae-c897edd4e447
md"""
#### With interference
"""

# ╔═╡ c9418469-8030-4f2d-bbd6-ce0e532170cb
function A3b(σs, ν; p)
	j0 = 1
	BWs = map(σ->BW(σ, p.m, p.Γ), σs)
	sum(BWs[1]*kronecker(ν,λ)*A(cosθ23(σs, ms^2),p.s,λ) + 
		BWs[2]*wignerd(j0,ν,λ, cosζ12_for0(σs, ms^2))*A(cosθ31(σs, ms^2),p.s,λ) *
			phase(2ν,2λ) +
		BWs[3]*wignerd(j0,ν,λ, cosζ31_for0(σs, ms^2))*A(cosθ12(σs, ms^2),p.s,λ)
			for λ in -j0:j0)
end

# ╔═╡ ae5d5bcf-6b51-41df-8b67-cd03611c509f
I3b(σs; p) = sum(abs2, A3b(σs, ν; p) for ν in -1:1)

# ╔═╡ aa65878e-b76e-4c44-8d89-5bd9a99bc2b4
begin
	plot(layout=grid(1,2), size=(650,300), grid=false, frame=:box)
	# 	
	plot!(sp=1, σs->I3b(σs; p=(s=1, m=mR, Γ=ΓR)), ms, c=:viridis)
	plot!(sp=1, border12(ms), l=(2,:black), lab="", title="P-wave")
	# 	
	plot!(sp=2, σs->I3b(σs; p=(s=3, m=mR, Γ=ΓR)), ms, c=:viridis)
	plot!(sp=2, border12(ms), l=(2,:black), lab="", title="F-wave")
end

# ╔═╡ Cell order:
# ╟─8dba788b-f74f-441c-a712-6116761a4a83
# ╠═3d7df13c-75a7-46ac-b28e-18f263de6af7
# ╠═66d3aa29-51c8-49bd-8df3-bc6469d18901
# ╠═54f6f89e-f496-4ad9-8047-44e9dbd43ef4
# ╠═a8c0f283-92fe-4433-b5dc-341b14c7a0a2
# ╟─9dcf87e4-beff-4b66-9118-6b6656e35e46
# ╟─876e7087-3948-4f84-9dcd-a4156c563983
# ╠═86ba31c9-57bf-4cd2-8eb0-d66be0ddc89f
# ╠═e3206b0d-277a-44e8-a26b-355aeb1dcae2
# ╠═b60ec657-a8bd-4229-9848-4f3a7fab7fde
# ╠═fc53f25c-c6d7-402c-9d1f-2f41daed0fe1
# ╟─69a05a98-f0e6-49de-be57-42b51dd0b551
# ╟─baa08139-93f6-4c6f-9fae-c897edd4e447
# ╠═c9418469-8030-4f2d-bbd6-ce0e532170cb
# ╠═ae5d5bcf-6b51-41df-8b67-cd03611c509f
# ╠═aa65878e-b76e-4c44-8d89-5bd9a99bc2b4
