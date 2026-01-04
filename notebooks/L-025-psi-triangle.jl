### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 342b5e0e-0ede-11f0-33ec-89efc7931f8a
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(url = "https://github.com/mmikhasenko/TriangleSingularity.jl.git"),
        Pkg.PackageSpec("Plots"),
        Pkg.PackageSpec("Parameters"),
        Pkg.PackageSpec("DelimitedFiles"),
        Pkg.PackageSpec("ThreeBodyDecays"),
        Pkg.PackageSpec("HadronicLineshapes"),
    ])
    # 
    using TriangleSingularity
    using HadronicLineshapes
    using ThreeBodyDecays
    using DelimitedFiles
    using Parameters
    using Plots
end

# ╔═╡ ee339c44-ac70-4315-8504-d37402d073dc
md"""
# Triangle mechanism for $T_{cc1}(4430)$

(M.Mikhasenko, RUB, 1.04.2025)

Triangle mechnism implies probability migration from one final state from the other.
For the case we consider, the re-distribution of probability is assumed to happen between the two final states:

```math
\begin{align}
B^0 &\to \psi\,\pi^+ K^-\\
B^0 &\to \psi(4260)\,\pi^+ K^-
\end{align}
```

The model also assumes:
1. Decay rate of $B^0 \to \psi(4260)\,\pi^+ K^-$ is significant (much larger that $B^0 \to \psi\,\pi^+ K^-$)
2. Coupled-channel between $\psi\,\pi^+$ and $\psi(4260)\,\pi^+$ is significant

### References for the triangle calculations:

- A book by Gribov [inspire:2180771](https://inspirehep.net/literature/2180771)
- Nature of a1(1420) [inspire:1341619](https://inspirehep.net/literature/1341619)
- Triangle Pentaquarks [inspire:1384521](https://inspirehep.net/literature/1384521)
- Narrow Pc (LHCb) [inspire:1728691](https://inspirehep.net/literature/1728691)

"""

# ╔═╡ a77fe681-4037-424c-ad21-d32e82f88b67
theme(:boxed, fontfamily = "Computer Modern")

# ╔═╡ 090b83c5-b302-45a8-8d3f-7acf2cd79d37
begin
    const m_D2600 = 2.62
    const Γ_D2600 = 0.14

    const mD = 1.8648
    const mπ = 0.13957
    const mψ = 3.686
    const m_ψ4260 = 4.222
    const Γ_ψ4260 = 0.049
    const mthr = m_ψ4260 + mπ

    const mKx = 0.8955
    const ΓKx = 0.047

    const mB = 5.279
    const mK = 0.494
end;

# ╔═╡ 052ed5c9-783a-48cf-8bc3-1f258de39633
const support = (mψ + mπ, 5.5)

# ╔═╡ e7e5a97c-0328-4d8d-b7a4-e70dddf41fa2
md"""
```
              o --- M₂
          m₁ / \
            /   \ m₃
    M₃ --- o --- o --- M₁
              m₂
```
"""

# ╔═╡ 49068b63-6c6d-4804-9a23-d2afcf95ff1a
function A_triangle(s)
    m1² = (mKx - 1im * ΓKx / 2)^2# 
    m2² = (m_ψ4260 - 1im * Γ_ψ4260 / 2)^2 # 
    m3² = mπ^2
    M3² = mB^2
    M2² = mK^2
    M1² = s
    return triangleloop(m1², m2², m3², M1², M2², M3²)
end

# ╔═╡ 821a6ecf-99cf-4135-877e-a71fb20af9f4
BW_Z = BreitWigner(4.41, 0.25)

# ╔═╡ ceefe3f4-120a-48e6-8d0a-14467a1b6692
const scale = A_triangle(BW_Z.m^2) / BW_Z(BW_Z.m^2)

# ╔═╡ 684e722b-91f0-46f5-a09a-45dbaa6989a9
BW_Z_scaled = BW_Z * (s -> scale);

# ╔═╡ ee227703-3671-4430-9981-6bc1c366d811
let
    plot(title = "Triangle versus Breit-Wigner")
    plot!(m -> imag(A_triangle(m^2)), support..., lw = 2, lab = "Triangle")
    plot!(m -> real(A_triangle(m^2)), support..., lw = 2)
    vspan!([support[1], m_ψ4260 + mπ], alpha = 0.2)
    #
    plot!(m -> imag(BW_Z_scaled(m^2)), support..., c = 1, alpha = 0.3, lab = "BW")
    plot!(m -> real(BW_Z_scaled(m^2)), support..., c = 2, alpha = 0.3)
    plot!(xlab = "\$m(\\psi\\pi)\$ [GeV]", ylab = "\$\\Re\\,\\mathcal{A}\$,  \$\\Im\\,\\mathcal{A}\$")
end

# ╔═╡ 1725193d-7262-4fa6-b9c4-4467c141053e
let
    plot(title = "Triangle versus Breit-Wigner")
    plot(m -> abs2(A_triangle(m^2)), support..., lw = 2, lab = "triangle")
    plot!(m -> abs2(BW_Z_scaled(m^2)), support..., c = 1, alpha = 0.3, lab = "BW")
    plot!(xlab = "\$m(\\psi\\pi)\$ [GeV]", ylab = "\$|\\mathcal{A}|^2\$")
end

# ╔═╡ d2fe0ef3-527a-477b-bc75-ccd1ea467c90
let
    mv = range(support..., 3000)
    Av = @. A_triangle(mv^2)
    plot(aspect_ratio = 1, leg = :topleft)
    plot!(Av, arrow = true, lw = 2, lab = "triangle")
    scatter!([A_triangle(mthr^2)], m = (4, :red, :o), lab = "nominal threshold")
    # 
    Bv = map(BW_Z_scaled, mv .^ 2)
    plot!(Bv, arrow = true, c = 1, alpha = 0.3, lab = "BW")

    extr = extrema(real.(vcat(Av, Bv)))
    xlim = extr .+ [-1, 1] * 0.05 * diff(collect(extr))[1]
    plot!(; xlim, xlab = "\$\\Re\\,\\mathcal{A}\$", ylab = "\$\\Im\\,\\mathcal{A}\$")
end

# ╔═╡ 0cac0956-db52-45aa-827a-52b18ab5e94d
let
    mv = range(support..., 100)
    Av = A_triangle.(mv .^ 2)
    data = [mv real.(Av) imag.(Av)]
end

# ╔═╡ c64817d1-cf89-4b6b-8218-36adf4a09400
md"""
## Kinematics

Two reactions overlap in the phase space allowing probability to leak from one final state to the other.
"""

# ╔═╡ 3869eb43-53ab-414d-8774-4cd6894d4a3e
ms = ThreeBodyMasses(mψ, mπ, mK; m0 = mB)

# ╔═╡ 6aa344c0-6411-4842-a91f-951a8ac14d57
ms_coupled = ThreeBodyMasses(m_ψ4260, mπ, mK; m0 = mB)

# ╔═╡ 5e383b7d-ecc2-424a-aec6-c12ae97f376d
begin
    plot(title = "\$B^0\$ decay phase space", aspect_ratio = 2)
    plot!(border31(ms) |> Shape, fillcolor = 3, fillalpha = 0.1)
    plot!(border31(ms_coupled) |> Shape, fillcolor = 2, fillalpha = 0.7)
    hline!([mKx^2], c = 4, lw = 5, ann = (4^2, mKx^2, text("\$K^\\star\$", :bottom)))
    plot!(xlab = "\$m^2(\\psi π)\$ [GeV\$^2\$]", ylab = "\$m^2(Kπ)\$ [GeV\$^2\$]")
end

# ╔═╡ Cell order:
# ╟─ee339c44-ac70-4315-8504-d37402d073dc
# ╟─342b5e0e-0ede-11f0-33ec-89efc7931f8a
# ╠═a77fe681-4037-424c-ad21-d32e82f88b67
# ╠═090b83c5-b302-45a8-8d3f-7acf2cd79d37
# ╠═052ed5c9-783a-48cf-8bc3-1f258de39633
# ╟─e7e5a97c-0328-4d8d-b7a4-e70dddf41fa2
# ╠═49068b63-6c6d-4804-9a23-d2afcf95ff1a
# ╠═821a6ecf-99cf-4135-877e-a71fb20af9f4
# ╠═ceefe3f4-120a-48e6-8d0a-14467a1b6692
# ╠═684e722b-91f0-46f5-a09a-45dbaa6989a9
# ╟─ee227703-3671-4430-9981-6bc1c366d811
# ╟─1725193d-7262-4fa6-b9c4-4467c141053e
# ╟─d2fe0ef3-527a-477b-bc75-ccd1ea467c90
# ╠═0cac0956-db52-45aa-827a-52b18ab5e94d
# ╟─c64817d1-cf89-4b6b-8218-36adf4a09400
# ╠═3869eb43-53ab-414d-8774-4cd6894d4a3e
# ╠═6aa344c0-6411-4842-a91f-951a8ac14d57
# ╟─5e383b7d-ecc2-424a-aec6-c12ae97f376d
