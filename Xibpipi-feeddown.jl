### A Pluto.jl notebook ###
# v0.19.12

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

# ╔═╡ adad5030-bce6-11ed-3d9e-e96ce1b40926
# ╠═╡ show_logs = false
begin
    # cd(mktempdir())
    using Pkg
    Pkg.add([PackageSpec("Plots"),
        PackageSpec("PlutoUI"),
        PackageSpec(
            url="https://github.com/mmikhasenko/ThreeBodyDecay.jl")])
    # 
    using Plots
    using ThreeBodyDecay
    using PlutoUI
end

# ╔═╡ b3dea2b0-8240-45cf-ab89-6f174a1767d0
theme(:wong2, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing, lab="",
    legendfontsize=9, legend=:topright)

# ╔═╡ a14ec568-7b0d-458e-9c51-fdc9d4ecc0a5
begin
    const mΞb⁰ = 5791.9
    const mΞbᶜ = 5797.0
    # 
    const mπ⁰ = 134.9768
    const mπᶜ = 139.57039
    # 
    const mΞb′ᶜ = 5935.13
    # 
    const mΞbˣᶜ = 5955.74
    const mΞbˣ⁰ = 5952.37
    # 
    const mΞbˣˣ1h⁰ = 6087.24
    # 
    const mΞbˣˣ3hᶜ = 6099.74
    const mΞbˣˣ3h⁰ = 6095.36
end

# ╔═╡ bc5932c1-4464-4df8-b80d-be042860b371
begin # isospin breaking
    ΔmΞb = mΞbᶜ - mΞb⁰
    ΔmΞbˣ = mΞbˣᶜ - mΞbˣ⁰
    ΔmΞbˣˣ3h = mΞbˣˣ3hᶜ - mΞbˣˣ3h⁰
    (; ΔmΞb, ΔmΞbˣ, ΔmΞbˣˣ3h)
end

# ╔═╡ 8f3fb833-246e-4905-a94b-ec4cb1cac8a7
md"""
## Decay paths
"""

# ╔═╡ db4565ba-cf0d-410e-8991-ce26418be5be
begin
    struct Particle
        symbol::Symbol
        mass::Float64
        charge::Int
    end
    Base.show(io::IO, p::Particle) = show(io, p.symbol)
    # 
    charge(p::Particle) = p.charge
    charge(c::Tuple) = sum(charge, c) |> isodd |> Int
    # 
    mass(p::Particle) = p.mass
    mass(c::Tuple) = sum(mass, c)
end

# ╔═╡ 1765da46-bcc2-4985-9313-4cf967be7092
function isrealistic((a, b, c))
    (charge(a) == charge(b) == charge(c)) &&
        (mass(a) > mass(b) > mass(c)) &&
    (b[2] == c[2] || b[2] == c[3])
end

# ╔═╡ 6353d970-bcf0-4749-8294-33abc38e2f4d
md"""
## Three-body spectra
"""

# ╔═╡ 5110f09a-6c39-4f91-bd8c-371c92c1410c
md"""
### with :Ξbˣ⁰

The spectrum has only one state, 3/2-
"""

# ╔═╡ d69361a5-ae3c-4ed8-adb6-2f82f0f8a238
md"""
### With Ξb′ᶜ

The spectrum can have two peaks, both 1/2, and 3/2
"""

# ╔═╡ 91d297c7-2a22-4aa9-80b9-0068cd19a80e
md"""
### with Ξbˣᶜ

The spectrum has only one peak, 3/2-
"""

# ╔═╡ cc6d1381-a039-45f1-a31a-908d3a057ed0
md"""
## Feeddown
"""

# ╔═╡ c488a4f4-525a-4504-a23d-f6cfe06ce291
md"""
## Compute the feeddown range
"""

# ╔═╡ 34813e68-05d1-4e9c-9843-f16725e90a2a
"""
	feeddown13((a,b,c))

Computes the mass of the (c1,c3) interval when the resonance is in (c1,c2)
"""
function feeddown13((a, b, c))
    mR = mass(b[1]) # resonance in [c1,c2]
    ms = ThreeBodyMasses(mass.(c)..., m0=(mass(a)))
    msq = ms^2
    σ2minmax = σ2of3.([-1, 1], Ref(mR^2), Ref(msq))
    return sqrt.(σ2minmax)
end

# ╔═╡ 972f828a-a08c-462f-bb87-f84bf6bd1b01
feeddownQ13(path) = feeddown13(path) .- mass((path[3][[1, 3]]...,))

# ╔═╡ a992d625-c3a1-4604-9686-8eeb8173b3c5
function isSwave((a,b,c))
	isspinhalf_Ξbˣˣ = string(a[1].symbol)[8]=='1'
	isspinhalf = string(b[1].symbol)[4]=='′'
	return isspinhalf_Ξbˣˣ == isspinhalf
end

# ╔═╡ 7756e2ea-05b0-478a-935d-b68bae55a3f7
md"""
### (Ξbᶜ,πᶜ) spectrum
"""

# ╔═╡ 1a6348ae-5639-41ed-ac95-75f0d44c19ef
md"""
### (Ξb⁰,πᶜ) spectrum
"""

# ╔═╡ 9aae6f45-08dd-4f9c-b943-0d9174c7db31
md"""
#### Masses of unobserved states

Mass difference $\Xi_b^{\prime-} - \Xi_b^{\prime0}$ : $(@bind shift Slider(3 : 0.01 : 6;	default=ΔmΞbˣ, show_value=true)) MeV

Mass difference $\Xi_b^{**-}(1/2^-) - \Xi_b^{**0}(1/2^-)$ :
	$(@bind shiftˣˣ Slider(3 : 0.01 : 6;
	default=(ΔmΞb + ΔmΞbˣ + ΔmΞbˣˣ3h)/3, show_value=true)) MeV

"""

# ╔═╡ 5c1752a2-5e99-4841-975c-1a2b19ad8fd7
begin
    mΞb′⁰ = mΞb′ᶜ - shift
    mΞbˣˣ1hᶜ = mΞbˣˣ1h⁰ + shiftˣˣ
end;

# ╔═╡ 9fb28a38-c3be-41a4-b7f8-1ade01df51a6
begin
    plot(layout=grid(1, 2), size=(500, 700), link=:y, title=["neutral" "charged"])
    # 
    [mΞb⁰ mΞb′⁰ mΞbˣ⁰ mΞbˣˣ1h⁰ mΞbˣˣ3h⁰] |>
    x -> plot!(sp=1, complex.(vcat((x .* 1im .+ [-1, 1, NaN])...)), lw=2)
    # 
    [mΞbᶜ mΞb′ᶜ mΞbˣᶜ mΞbˣˣ1hᶜ mΞbˣˣ3hᶜ] |>
    x -> plot!(sp=2, complex.(vcat((x .* 1im .+ [-1, 1, NaN])...)), lw=2)
    # 
    # 
    [mΞb⁰ + mπ⁰ mΞbᶜ + mπᶜ] |> # mΞb⁰+2mπ⁰  mΞb⁰+2mπᶜ] |>
    x -> plot!(sp=1, complex.(vcat((x .* 1im .+ [-1, 1, NaN])...)), ls=:dash)
    [mΞb′⁰ + mπ⁰ mΞb′ᶜ + mπᶜ] |>
    x -> plot!(sp=1, complex.(vcat((x .* 1im .+ [-1, 1, NaN])...)), ls=:dashdot)
    [mΞbˣ⁰ + mπ⁰ mΞbˣᶜ + mπᶜ] |>
    x -> plot!(sp=1, complex.(vcat((x .* 1im .+ [-1, 1, NaN])...)), ls=:dashdot)
    # 
    # 
    [mΞbᶜ + mπ⁰ mΞb⁰ + mπᶜ] |> # mΞbᶜ+2mπ⁰  mΞbᶜ+2mπᶜ] |>
    x -> plot!(sp=2, complex.(vcat((x .* 1im .+ [-1, 1, NaN])...)), ls=:dash)
    [mΞb′ᶜ + mπ⁰ mΞb′⁰ + mπᶜ] |>
    x -> plot!(sp=2, complex.(vcat((x .* 1im .+ [-1, 1, NaN])...)), ls=:dashdot)
    [mΞbˣᶜ + mπ⁰ mΞbˣ⁰ + mπᶜ] |>
    x -> plot!(sp=2, complex.(vcat((x .* 1im .+ [-1, 1, NaN])...)), ls=:dashdot)
    # 
    plot!(xaxis="", yaxis="")
end

# ╔═╡ 30ed523b-c969-4c8b-bdbd-e5f65ce0f636
begin
    Ξbˣˣ1hᶜ = Particle(:Ξbˣˣ1hᶜ, mΞbˣˣ1hᶜ, -1)
    Ξbˣˣ1h⁰ = Particle(:Ξbˣˣ1h⁰, mΞbˣˣ1h⁰, 0)
    Ξbˣˣ3hᶜ = Particle(:Ξbˣˣ3hᶜ, mΞbˣˣ3hᶜ, -1)
    Ξbˣˣ3h⁰ = Particle(:Ξbˣˣ3h⁰, mΞbˣˣ3h⁰, 0)
    # 
    Ξb′⁰ = Particle(:Ξb′⁰, mΞb′⁰, 0)
    Ξb′ᶜ = Particle(:Ξb′ᶜ, mΞb′ᶜ, -1)
    Ξbˣ⁰ = Particle(:Ξbˣ⁰, mΞbˣ⁰, 0)
    Ξbˣᶜ = Particle(:Ξbˣᶜ, mΞbˣᶜ, -1)
    #
    Ξb⁰ = Particle(:Ξb⁰, mΞb⁰, 0)
    Ξbᶜ = Particle(:Ξbᶜ, mΞbᶜ, -1)
    # 
    π⁰ = Particle(:π⁰, mπ⁰, 0)
    πᶜ = Particle(:πᶜ, mπᶜ, -1)
end ;

# ╔═╡ 97cda745-dbd8-4138-8767-7806bb7a8303
allpaths = let
	xx_states = ((Ξbˣˣ1hᶜ,), (Ξbˣˣ3hᶜ,), (Ξbˣˣ1h⁰,), (Ξbˣˣ3h⁰,))
	xpi_states = [(x, y) for x in [Ξb′⁰, Ξb′ᶜ, Ξbˣ⁰, Ξbˣᶜ]
	              for y in (π⁰, πᶜ)]
	pipi_states = let pi = (π⁰, πᶜ)
	    [(x, pi[y], pi[z]) for x in [Ξb⁰, Ξbᶜ]
	     for y in 1:2 for z in y:2]
	end
	# combine 
	[(xx, xpi, pipi) for xx in xx_states
	            for xpi in xpi_states for pipi in pipi_states]
end ;

# ╔═╡ 2d4b0498-d05e-4927-a21f-2eff73fce3f2
realisticpaths = filter(isrealistic, allpaths) ;

# ╔═╡ c04ade76-3a10-4527-852c-adff841424c2
measurablepaths = filter(realisticpaths) do (a, b, c)
    !(π⁰ ∈ c[2:3])
end

# ╔═╡ f3ef454a-0e54-4f99-8060-ebf060fe3821
filter(measurablepaths) do (a, b, c)
    b[1] == Ξbˣ⁰
end

# ╔═╡ 8a60b609-b4ba-4dfa-847d-59b8b7e79f49
filter(measurablepaths) do (a, b, c)
    b[1] == Ξb′ᶜ
end

# ╔═╡ 67221e2c-a389-4db3-87f9-c8e98acfb3fc
filter(measurablepaths) do (a, b, c)
    b[1] == Ξbˣᶜ
end

# ╔═╡ 32d7821c-b3bc-4af6-a3bd-93a88a5829c2
f1 = filter(realisticpaths) do (a, b, c)
    c[1] == Ξbᶜ &&
        πᶜ ∈ c[2:3] &&
        !(b[1] ∈ (Ξbˣ⁰, Ξb′⁰))
end

# ╔═╡ 79a88a33-bd5c-464e-a237-e82471fba555
let
    plot(xlim=(0, 50), ylim=(0, 1.5))
    # Qranges = feeddownQ13.(f1)
    [Ξb′⁰, Ξbˣ⁰] .|>
    m -> plot!((mass(m) - mΞbᶜ - mπᶜ) + 1im .+ [-1, 1];
        fill=0, alpha=0.5, lab=string(m))
    # Qranges .|> Qs->plot!(Qs .+ 0.5im, fill=0, alpha=0.3)
    f1 .|> path -> plot!(feeddownQ13(path) .+ (1im*(1+isSwave(path))/4), fill=0, alpha=0.3,
        lab=string(path[1][1]) * "→ $(path[2][1])(→$(path[3][1]) $(path[3][2])) $(path[3][3])")
    plot!(xlab="m(Ξbᶜ πᶜ)-m(Ξbᶜ)-m(πᶜ) (MeV)", ylab="")
end

# ╔═╡ 77a5446d-818a-4000-8ca6-5dd7f3c08b80
f2 = filter(realisticpaths) do (a, b, c)
    c[1] == Ξb⁰ &&
        (c[2] == πᶜ || c[3] == πᶜ) &&
        !(b[1] ∈ (Ξbˣᶜ, Ξb′ᶜ))
end

# ╔═╡ 25bddd7f-61d8-4a77-baef-98ca620ef818
let
    plot(xlim=(0, 50), ylim=(0, 1.5))
    [Ξb′ᶜ, Ξbˣᶜ] .|>
    m -> plot!((mass(m) - mΞb⁰ - mπᶜ) + 1im .+ [-1, 1];
        fill=0, alpha=0.5, lab=string(m))
    f2 .|> path -> plot!(feeddownQ13(path) .+ (1im*(1+isSwave(path))/4), fill=0, alpha=0.3,
        lab=string(path[1][1]) * "→ $(path[2][1])(→$(path[3][1]) $(path[3][2])) $(path[3][3])")
    plot!(xlab="m(Ξb⁰ πᶜ)-m(Ξb⁰)-m(πᶜ) (MeV)", ylab="")
end

# ╔═╡ Cell order:
# ╠═adad5030-bce6-11ed-3d9e-e96ce1b40926
# ╠═b3dea2b0-8240-45cf-ab89-6f174a1767d0
# ╠═a14ec568-7b0d-458e-9c51-fdc9d4ecc0a5
# ╠═bc5932c1-4464-4df8-b80d-be042860b371
# ╠═5c1752a2-5e99-4841-975c-1a2b19ad8fd7
# ╠═9fb28a38-c3be-41a4-b7f8-1ade01df51a6
# ╟─8f3fb833-246e-4905-a94b-ec4cb1cac8a7
# ╠═db4565ba-cf0d-410e-8991-ce26418be5be
# ╠═30ed523b-c969-4c8b-bdbd-e5f65ce0f636
# ╠═97cda745-dbd8-4138-8767-7806bb7a8303
# ╠═1765da46-bcc2-4985-9313-4cf967be7092
# ╠═2d4b0498-d05e-4927-a21f-2eff73fce3f2
# ╠═c04ade76-3a10-4527-852c-adff841424c2
# ╟─6353d970-bcf0-4749-8294-33abc38e2f4d
# ╟─5110f09a-6c39-4f91-bd8c-371c92c1410c
# ╠═f3ef454a-0e54-4f99-8060-ebf060fe3821
# ╟─d69361a5-ae3c-4ed8-adb6-2f82f0f8a238
# ╠═8a60b609-b4ba-4dfa-847d-59b8b7e79f49
# ╟─91d297c7-2a22-4aa9-80b9-0068cd19a80e
# ╠═67221e2c-a389-4db3-87f9-c8e98acfb3fc
# ╟─cc6d1381-a039-45f1-a31a-908d3a057ed0
# ╟─c488a4f4-525a-4504-a23d-f6cfe06ce291
# ╠═34813e68-05d1-4e9c-9843-f16725e90a2a
# ╠═972f828a-a08c-462f-bb87-f84bf6bd1b01
# ╠═a992d625-c3a1-4604-9686-8eeb8173b3c5
# ╟─7756e2ea-05b0-478a-935d-b68bae55a3f7
# ╠═32d7821c-b3bc-4af6-a3bd-93a88a5829c2
# ╠═79a88a33-bd5c-464e-a237-e82471fba555
# ╟─1a6348ae-5639-41ed-ac95-75f0d44c19ef
# ╠═77a5446d-818a-4000-8ca6-5dd7f3c08b80
# ╟─25bddd7f-61d8-4a77-baef-98ca620ef818
# ╟─9aae6f45-08dd-4f9c-b943-0d9174c7db31
