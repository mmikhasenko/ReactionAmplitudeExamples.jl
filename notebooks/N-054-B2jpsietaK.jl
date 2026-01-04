### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 545b55d0-d235-11ee-25ce-8d43a2e4befa
# ╠═╡ show_logs = false
begin
    using Pkg: add, PackageSpec, activate
    activate(mktempdir())
    # 
    # does not work
    # PackageSpec(; path = "/Users/mikhailmikhasenko/Documents/PDG.CAT/PDGdb.jl")
    add([
        PackageSpec(; url="https://github.com/mmikhasenko/PDGdb.jl.git"),
        PackageSpec("Downloads"),
        PackageSpec("Plots"),
        PackageSpec("Parameters"),
        PackageSpec("ThreeBodyDecays"),
        PackageSpec("Measurements"),
        PackageSpec("LaTeXStrings")])
    # 
    using ThreeBodyDecays
    using Parameters
    using PDGdb
    using PDGdb.DataFrames
    using Downloads
    using Plots
    using Plots.PlotMeasures: mm
    using Measurements
    using LaTeXStrings
end

# ╔═╡ 87c869d6-74fd-42b2-a3d6-819242e3eeb7
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto),
    lw=1.2, lab="", colorbar=false)

# ╔═╡ 486bb34f-6911-4ba8-9bd0-e20e69f0eed2
# m0->0,0,0
# σ1 ∈ [0, m0]ˆ2
# σ3 ∈ [0, m0]ˆ2

# Phi3 = 1/s ∫ dσ1 dσ3 = 1/s * [m0ˆ2 * m0ˆ2/2]->m0ˆ2

# phase space: m0ˆ2

# m^2*(m/2-E)*E

# ╔═╡ 79cd619e-c28d-4d66-a0f9-9204dbabb78f


# ╔═╡ c1af5512-c2a3-4fda-aef9-de777cd0bf1b
# ╠═╡ show_logs = false
begin
    db_url = "https://pdg.lbl.gov/2023/api/pdg-2023-v0.0.5.sqlite"
    tmp_dir = mktempdir()
    const db_path = joinpath(tmp_dir, "pdg-2023-v0.0.5.sqlite")
    Downloads.download(db_url, joinpath(tmp_dir, "pdg-2023-v0.0.5.sqlite"))
    PDGdb.connect(db_path)
end;

# ╔═╡ d2a87e31-e98b-4842-bc3f-fbff84e907d4
begin
    mηp = (PDGdb.pdg("eta^'(958)")|>properties|>mass).value[1] ./ 1e3
    mJψ = ((PDGdb.pdg("J/psi(1S)")|>properties|>mass).value[1] ./ 1e3)
    mK = (PDGdb.pdg("K+")|>properties|>mass).value[1] ./ 1e3
    mB = (PDGdb.pdg("B+")|>properties|>mass).value[1] ./ 1e3
end

# ╔═╡ cf7d38a3-a27a-40fa-b73e-327e51df5aff
ms = ThreeBodyMasses(mJψ, mηp, mK; m0=mB);

# ╔═╡ 39c21ef6-297c-481f-9b88-e74f0f96eb1d
tbs = ThreeBodySystem(ms, ThreeBodySpins(2, 0, 0; two_h0=0))

# ╔═╡ bd8abd3c-fbba-4d92-adc9-496f4e899e65
mΓ1430 = NamedTuple{(:m, :Γ)}(
    (PDGdb.pdg("K_0^*(1430)0") |> properties |> parameters).value ./ 1e3)

# ╔═╡ 1f354ff4-adcb-4d36-9935-351ec78806a1
begin
    @with_kw struct BW
        m::Float64
        Γ::Float64
    end
    (bw::BW)(σ) = 1 / (bw.m^2 - σ - 1im * bw.m * bw.Γ)
end

# ╔═╡ 62ab3830-e8f7-4878-bae1-aa54d4b96245
PC, PV = ['-', '-', '-', '-'], ['-', '-', '-', '+']

# ╔═╡ dabf060b-dc6c-47d5-bc72-b84da7e87b64
model = let
    chains = map([PC, PV]) do Ps
        DecayChainsLS(;
            k=1, Xlineshape=BW(; mΓ1430...), two_j=2, parity='-',
            Ps, tbs) |> vec
    end |> x -> vcat(x...)
    d = ThreeBodyDecay("K_0^*(1430)0" .=> zip([-3.0im, 4.0, 0.0], chains))
end;

# ╔═╡ 4cec4ef7-1f0d-4826-8b7a-16ceeeeb9245
sqrt_nt(x) = typeof(x)(sqrt.(Tuple(x)))

# ╔═╡ 7979f78c-3cfc-4feb-bfd3-032f38ddf323
begin
    heatmap(range(4, 4.8, 100), range(1.4, 2.25, 100), (e3, e1) -> begin
            σs = Invariants(ms; σ1=e1^2, σ3=e3^2)
            !inphrange(σs, ms) && return NaN
            unpolarized_intensity(model, σs) * e1 * e3
        end, c=cgrad(:lajolla, rev=true))
    # 
    plot!(sqrt_nt.(border31(ms)), lw=3, lc=:black,
        ylab="m²(K⁺ η') [GeV]",
        xlab="m²(J/ψ η') [GeV]",
        size=(500, 400))
end

# ╔═╡ 129d612d-f347-4efe-9853-91aeb8b58690
border31(ms)

# ╔═╡ b4bd3722-e7fb-41df-840b-6314a04fed8e
mth =
    (PDGdb.pdg("Lambda_c()+")|>properties|>mass).value[1] ./ 1e3 +
    (PDGdb.pdg("D0")|>properties|>mass).value[1] ./ 1e3

# ╔═╡ ed32267e-18c0-4130-8fda-47b1b39100db
[4.3, 4.5] .- mth

# ╔═╡ 7dee02e2-a418-47bf-938a-e33ad9967780
mthX = (PDGdb.pdg("Lambda_c()+")|>properties|>mass).value[1] ./ 1e3 +
       (PDGdb.pdg("D^*()-")|>properties|>mass).value[1] ./ 1e3

# ╔═╡ 2d968856-b561-46f7-88ad-21efe970bd32
[4.3, 4.5] .- mthX

# ╔═╡ 13428f1c-af0e-4df6-a398-92c9d76c99e9
(PDGdb.pdg("Lambda_c()+")|>properties|>mass).value[1] ./ 1e3 +
(PDGdb.pdg("D0")|>properties|>mass).value[1] ./ 1e3 +
(PDGdb.pdg("pi+")|>properties|>mass).value[1] ./ 1e3

# ╔═╡ b9ca5156-7f14-4836-ab54-016e04c7421b
(4.14 ± 0.45) / (4.3 ± 1.3)

# ╔═╡ 9b1e620a-0806-4ade-9a9b-fc6dc65f04af
mc = 1.5

# ╔═╡ 394e032d-2430-4443-97ea-f7f87f8f9608
mb = 4.0

# ╔═╡ 76a5ff6f-6087-47b2-80a1-256d9d5f6759
ρ(m, m1, m2) = sqrt(m^2 - (m1 + m2)^2) * sqrt(m^2 - (m1 + m2)^2) / m^2

# ╔═╡ f5cb6263-cc05-4f3f-b5e1-17dd0c75c668
begin
    mBc = (PDGdb.pdg("B_c()+")|>properties|>mass).value[1] ./ 1e3
    mBs = (PDGdb.pdg("B_s()0")|>properties|>mass).value[1] ./ 1e3
    mπ = (PDGdb.pdg("pi+")|>properties|>mass).value[1] ./ 1e3
end

# ╔═╡ e328cf02-ac34-428b-a46e-89f495464ebe
V_bc = 0.04

# ╔═╡ 5f642921-3ee9-47e2-ab97-1008758066f5
Br_c = mc^5 * ρ(mBc, mJψ, mπ)

# ╔═╡ 4bc34d1a-01f6-42d4-bbe2-595705245c96
Br_b = V_bc^2 * mb^5 * ρ(mBc, mBs, mπ)

# ╔═╡ 9c3c9900-7143-46ef-be52-a450602194fa
(mc)^5 / (V_bc^2 * (mb)^5)

# ╔═╡ 390ab0f4-0997-41f4-a606-4e8558bf16ae
Br_c / Br_b

# ╔═╡ 0f91c253-62a6-4bc0-8d79-b0f0989488e7
begin
    Bc2Bsπ = [167, 15.8, 58.4, 34.8]
    Bc2Jψπ = [1.43, 1.22, 1.97, 0.82]
    labels_Bc = ["Kiselev et al.", "Abd El-Hady et al.", "Chang et al.", "Anisimov et al."]
end

# ╔═╡ fe213793-269d-4559-bba6-52e66ae8cfcf
Bsπ_over_Jψπ = Bc2Bsπ ./ Bc2Jψπ

# ╔═╡ b70b8f8a-1ca8-4e35-a2b4-ec14dca020c8
16.4 / 0.13

# ╔═╡ 8be08df3-e19a-47ac-8909-48be38209716
Bsπ_over_Jψπ

# ╔═╡ e01bd2e7-f312-4c06-bbd7-e5433535410c
Bsπ_over_Jψπ_LHCb_syst_stat = 91 ± sqrt(10^2 + 8^2 + 3^2)

# ╔═╡ c92d1051-aea5-4a17-84aa-90e802ccf793
Bsπ_over_Jψπ_LHCb_stat = 91 ± 10

# ╔═╡ 032edc41-0a97-416a-83ab-8e5940bf05e9
begin
    plot(xlims=(-0.2, 4.3), ylim=(0, 130), size=(600, 150))
    scatter!(1:4, Bsπ_over_Jψπ, m=(6, :r), mc=4)
    hspan!(
        Bsπ_over_Jψπ_LHCb_syst_stat.val .+
        [-1, 1] .* Bsπ_over_Jψπ_LHCb_syst_stat.err, α=0.4, c=:orange)
    hspan!(Bsπ_over_Jψπ_LHCb_stat.val .+
           [-1, 1] .* Bsπ_over_Jψπ_LHCb_stat.err, α=0.2, c=:yellow)
    # 
    scatter!([0.0, 0.0], [Bsπ_over_Jψπ_LHCb_syst_stat, Bsπ_over_Jψπ_LHCb_stat], l=(2, :black), m=nothing)
    plot!(xticks=(0:4, ["LHCb", labels_Bc...]), ylab=L"\frac{\mathrm{B}(B_c^+\!\!\!\!\!\to B_s^0\pi^+\!\!)}{\mathrm{B}(B_c^+\!\!\!\!\!\to J\!/\!\psi\,\,\pi^+\!\!)}", left_margin=8mm)
end

# ╔═╡ e723a31e-ff14-46f4-83ca-648d694f0326
begin
    τDp = pick(PDGdb.pdg("D+") |> properties, [:pdgid => "S031T"])
    τDz = pick(PDGdb.pdg("D0") |> properties, [:pdgid => "S032T"])
    τBp = pick(PDGdb.pdg("B+") |> properties, [:pdgid => "S041T", :value_type => "AC"])
    τBz = pick(PDGdb.pdg("B0") |> properties, [:pdgid => "S042T", :value_type => "AC"])
    (; τDp, τDz, τBz, τBp)
end

# ╔═╡ 9e1db2e7-e63e-43f4-8833-f57708f34673
(τBz / τDz, τBp / τDp)

# ╔═╡ c18fadb1-e3f3-4eb8-8995-64415fa0e99f
1 / τDp

# ╔═╡ 5bc4418b-8fee-4651-9667-4d2555eb7ebf

# ╔═╡ Cell order:
# ╠═545b55d0-d235-11ee-25ce-8d43a2e4befa
# ╠═87c869d6-74fd-42b2-a3d6-819242e3eeb7
# ╠═486bb34f-6911-4ba8-9bd0-e20e69f0eed2
# ╠═79cd619e-c28d-4d66-a0f9-9204dbabb78f
# ╠═c1af5512-c2a3-4fda-aef9-de777cd0bf1b
# ╠═d2a87e31-e98b-4842-bc3f-fbff84e907d4
# ╠═cf7d38a3-a27a-40fa-b73e-327e51df5aff
# ╠═39c21ef6-297c-481f-9b88-e74f0f96eb1d
# ╠═bd8abd3c-fbba-4d92-adc9-496f4e899e65
# ╠═1f354ff4-adcb-4d36-9935-351ec78806a1
# ╠═62ab3830-e8f7-4878-bae1-aa54d4b96245
# ╠═dabf060b-dc6c-47d5-bc72-b84da7e87b64
# ╠═4cec4ef7-1f0d-4826-8b7a-16ceeeeb9245
# ╠═7979f78c-3cfc-4feb-bfd3-032f38ddf323
# ╠═129d612d-f347-4efe-9853-91aeb8b58690
# ╠═b4bd3722-e7fb-41df-840b-6314a04fed8e
# ╠═ed32267e-18c0-4130-8fda-47b1b39100db
# ╠═7dee02e2-a418-47bf-938a-e33ad9967780
# ╠═2d968856-b561-46f7-88ad-21efe970bd32
# ╠═13428f1c-af0e-4df6-a398-92c9d76c99e9
# ╠═b9ca5156-7f14-4836-ab54-016e04c7421b
# ╠═9b1e620a-0806-4ade-9a9b-fc6dc65f04af
# ╠═394e032d-2430-4443-97ea-f7f87f8f9608
# ╠═76a5ff6f-6087-47b2-80a1-256d9d5f6759
# ╠═f5cb6263-cc05-4f3f-b5e1-17dd0c75c668
# ╠═e328cf02-ac34-428b-a46e-89f495464ebe
# ╠═5f642921-3ee9-47e2-ab97-1008758066f5
# ╠═4bc34d1a-01f6-42d4-bbe2-595705245c96
# ╠═9c3c9900-7143-46ef-be52-a450602194fa
# ╠═390ab0f4-0997-41f4-a606-4e8558bf16ae
# ╠═0f91c253-62a6-4bc0-8d79-b0f0989488e7
# ╠═fe213793-269d-4559-bba6-52e66ae8cfcf
# ╠═032edc41-0a97-416a-83ab-8e5940bf05e9
# ╠═b70b8f8a-1ca8-4e35-a2b4-ec14dca020c8
# ╠═8be08df3-e19a-47ac-8909-48be38209716
# ╠═e01bd2e7-f312-4c06-bbd7-e5433535410c
# ╠═c92d1051-aea5-4a17-84aa-90e802ccf793
# ╠═e723a31e-ff14-46f4-83ca-648d694f0326
# ╠═9e1db2e7-e63e-43f4-8833-f57708f34673
# ╠═c18fadb1-e3f3-4eb8-8995-64415fa0e99f
# ╠═5bc4418b-8fee-4651-9667-4d2555eb7ebf
