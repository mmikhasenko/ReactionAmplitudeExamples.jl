### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ d8489d70-cf9f-11ec-1f26-0b616e422d52
# ╠═╡ show_logs = false
begin
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"),
        Pkg.PackageSpec("Plots"),
        Pkg.PackageSpec("LaTeXStrings"),
    ])
    # 
    using ThreeBodyDecay
    using Plots
    using LaTeXStrings
end

# ╔═╡ 47666396-b2dd-40c7-afe2-02fa25591963
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lw=1, lab="", colorbar=false,
    xlab=L"m^2_{\pi^+ \pi^-}\,\,(\mathrm{GeV}^2)",
    ylab=L"m^2_{\pi^+ \pi^-}\,\,(\mathrm{GeV}^2)")

# ╔═╡ e357f507-0730-4a4b-9c5e-689adf4734be
begin
    const mπ = 0.14
    const mρ, Γρ = 0.77, 0.15
    const ma1 = 1.26
end;

# ╔═╡ 21922eab-95bc-436e-8b1d-f532cdfa99e8
const ms = ThreeBodyMasses(mπ, mπ, mπ; m0=ma1)

# ╔═╡ fc6deed2-3dd4-4b82-9ea7-52d3ffd3af1e
const tbs = ThreeBodySystem(; ms, two_js=[0, 0, 0, 1] |> x2)

# ╔═╡ 19b9525b-6c0d-4677-9055-bc84e1dd3251
md"""
## Schematic Dalitz Plot
"""

# ╔═╡ 5c94caf5-7b2d-4142-a66e-abf1ad691989
begin
    plot(border12(ms), l=(1, :black), lab="")
    # 
    vspan!((mρ .+ Γρ / 2 .* [-1, 1]) .^ 2, α=0.3, c=2)
    vline!([mρ^2], lc=2, lw=3, ann=(0.3, mρ, text(L"\rho", 30)))
    # 
    hspan!((mρ .+ Γρ / 2 .* [-1, 1]) .^ 2, α=0.3, c=3)
    hline!([mρ^2], lc=3, lw=3, ann=(mρ, 0.3, text(L"\rho", 30)))
    # 
    plot!(border12(ms), l=(4, :black), lab="", aspectratio=1)
    #
    savefig("a1dalitzplotwithrho_schematic.pdf")
    plot!()
end

# ╔═╡ 94d75c16-0d2a-463b-b543-45da7ee08c24
md"""
## More realistic Dalitz Plot
"""

# ╔═╡ 260d3e58-e4bb-4512-b695-e8b682d19641
function BW(σ)
    p = sqrt(Kallen(σ, mπ^2, mπ^2) / (4σ))
    p0 = sqrt(Kallen(mρ^2, mπ^2, mπ^2) / (4σ))
    R = 1.5 # 1/GeV
    ff = sqrt(1 + R^2 * p0^2) / sqrt(1 + R^2 * p^2)
    Γ = Γρ * ff^2 * (p / p0)^3 * mρ / sqrt(σ)
    p * ff / (mρ^2 - σ - 1im * mρ * Γ)
end

# ╔═╡ e0e86f1e-c4d9-4558-80b7-42979f62ff7f
begin
    plot(σs -> abs2(BW(σs[1]) + BW(σs[2])), ms, aspectratio=1)
    plot!(border12(ms), l=(4, :black), lab="", aspectratio=1)
    # 
    savefig("a1dalitzplotwithrho_scalar.pdf")
    plot!()
end

# ╔═╡ c5ed4542-a60a-4d7e-bb57-ce9653ebfd6f
md"""
## Even-more realistic Dalitz Plot (+spin)
"""

# ╔═╡ 9926b8b2-e1a0-4cf6-b1b1-acd2c61dcd78
const dc = [DecayChain(; k, Xlineshape=BW, two_s=1 |> x2,
    Hij=RecouplingLS(; two_j=1 |> x2, two_ls=(1, 0) |> x2, two_ja=0, two_jb=0),
    HRk=RecouplingLS(; two_j=1 |> x2, two_ls=(0, 1) |> x2, two_ja=1 |> x2, two_jb=0),
    tbs) for k in 1:2]

# ╔═╡ 3e672231-c34c-4c99-8d68-2802607a20cc
A(x...) = amplitude(dc[1], x...) - amplitude(dc[2], x...)

# ╔═╡ 0cae4684-f9ef-4da2-8020-ec43e3b71078
const I = summed_over_polarization(abs2 ∘ A, tbs.two_js)

# ╔═╡ ab6185a0-f834-4023-b59d-0547e3fa01a5
begin
    plot(I, ms, aspectratio=1)
    plot!(border12(ms), l=(4, :black), lab="", aspectratio=1)
end

# ╔═╡ Cell order:
# ╠═d8489d70-cf9f-11ec-1f26-0b616e422d52
# ╠═47666396-b2dd-40c7-afe2-02fa25591963
# ╠═e357f507-0730-4a4b-9c5e-689adf4734be
# ╠═21922eab-95bc-436e-8b1d-f532cdfa99e8
# ╠═fc6deed2-3dd4-4b82-9ea7-52d3ffd3af1e
# ╟─19b9525b-6c0d-4677-9055-bc84e1dd3251
# ╠═5c94caf5-7b2d-4142-a66e-abf1ad691989
# ╟─94d75c16-0d2a-463b-b543-45da7ee08c24
# ╠═260d3e58-e4bb-4512-b695-e8b682d19641
# ╠═e0e86f1e-c4d9-4558-80b7-42979f62ff7f
# ╟─c5ed4542-a60a-4d7e-bb57-ce9653ebfd6f
# ╠═9926b8b2-e1a0-4cf6-b1b1-acd2c61dcd78
# ╠═3e672231-c34c-4c99-8d68-2802607a20cc
# ╠═0cae4684-f9ef-4da2-8020-ec43e3b71078
# ╠═ab6185a0-f834-4023-b59d-0547e3fa01a5
