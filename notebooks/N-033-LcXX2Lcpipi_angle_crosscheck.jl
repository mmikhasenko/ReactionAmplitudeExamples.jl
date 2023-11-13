### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 683556a0-c495-11ec-1724-85835878175b
begin
    import Pkg
    Pkg.activate(mktempdir())
    # 
    Pkg.add([
        Pkg.PackageSpec("Corpuscles"),
        Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl")
    ])
    using ThreeBodyDecay
    using Corpuscles
    using Corpuscles.Unitful
end

# ╔═╡ dfba203c-12cf-4e73-b25f-c7521e95b4e4
md"""
# Cross check of wigner angles
"""

# ╔═╡ 7735bcf3-a395-4ca3-bed9-3e9480fe437b
begin
    const mΛcˣˣ = convert(Float64,
        Particle("Lambda(c)(2625)").mass.value / 1u"GeV*c^-2")
    const mΛc = convert(Float64,
        Particle("Lambda(c)").mass.value / 1u"GeV*c^-2")
    const mπ = convert(Float64,
        Particle("pi").mass.value / 1u"GeV*c^-2")
end;

# ╔═╡ fce3386c-7e11-4690-8599-437d1ff5e280
const ms = ThreeBodyMasses(mΛc, mπ, mπ; m0=mΛcˣˣ)

# ╔═╡ 7c107c24-271c-447b-a15c-fb40544d3988
σs0 = randomPoint(ms)

# ╔═╡ 326e8fae-aa71-41d1-bd6e-46434ada1edf
allcosζ = (:cosζ21_for1, :cosζ21_for2, :cosζ13_for1, :cosζ13_for3, :cosζ32_for3, :cosζ32_for2, :cosζ12_for3, :cosζ23_for1, :cosζ31_for2, :cosζ12_for0, :cosζ23_for0, :cosζ31_for0)

# ╔═╡ b6d80036-2b4f-4bc2-b04a-4b0d982b308c
allcosθ = (:cosθ23, :cosθ31, :cosθ12)

# ╔═╡ d8526cb4-dfa3-41d4-91a8-0332cd2437c1
NamedTuple{allcosζ}([eval(:($(a)($(σs0), $(ms^2)))) for a in allcosζ])

# ╔═╡ c956662a-82b9-44a5-bb21-17f49182f58e
NamedTuple{allcosθ}([eval(:($(a)($(σs0), $(ms^2)))) for a in allcosθ])

# ╔═╡ Cell order:
# ╟─dfba203c-12cf-4e73-b25f-c7521e95b4e4
# ╠═683556a0-c495-11ec-1724-85835878175b
# ╠═7735bcf3-a395-4ca3-bed9-3e9480fe437b
# ╠═fce3386c-7e11-4690-8599-437d1ff5e280
# ╠═7c107c24-271c-447b-a15c-fb40544d3988
# ╠═326e8fae-aa71-41d1-bd6e-46434ada1edf
# ╠═b6d80036-2b4f-4bc2-b04a-4b0d982b308c
# ╠═d8526cb4-dfa3-41d4-91a8-0332cd2437c1
# ╠═c956662a-82b9-44a5-bb21-17f49182f58e
