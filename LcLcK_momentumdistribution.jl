### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 1125f9e0-d4ea-11ec-169b-d754f5700778
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.add([
		Pkg.PackageSpec(url="https://github.com/mmikhasenko/ThreeBodyDecay.jl"),
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("LaTeXStrings")
		])
	# 
	using ThreeBodyDecay
	using Plots
	using LaTeXStrings
	using QuadGK
end

# ╔═╡ 9e8a48f6-9e5b-4907-81cd-cabc1ac28213
md"""
# Check of $B^+ \to \Lambda_c^+ \overline{\Lambda}_c^-K^+$ efficiency
"""

# ╔═╡ 8b888bf9-7a8f-45e8-b268-41e9fdd8901e
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    lw=1.2, lab="", colorbar=false)

# ╔═╡ f9e829e0-6b08-403d-ba18-697ac51edfad
begin
	const mΛc = 2.28646
	const mK = 0.493677
	const mB = 5.27934
end ;

# ╔═╡ adfc0267-f6e6-412d-80db-3ee0c73bc32b
ms = ThreeBodyMasses(mΛc, mΛc, mK; m0=mB)

# ╔═╡ 0f28d98f-b7ec-4f98-abd7-017b5b5e4641
function project1(f, σ1, ms)
	function integrand(σ2)
		σs = Invariants(ms; σ1, σ2)
		return inphrange(σs,ms) ? f(σs) : 0.0
	end
	quadgk(integrand, lims2(ms)...)[1]
end

# ╔═╡ e2f48768-b995-4260-9e19-a1764af3082e
function project3(f, σ3, ms)
	function integrand(σ2)
		σs = Invariants(ms; σ2, σ3)
		return inphrange(σs,ms) ? f(σs) : 0.0
	end
	quadgk(integrand, lims2(ms)...)[1]
end

# ╔═╡ a802e131-2c78-4c3f-9d9f-49aa45f261ac
plot(map(x->sqrt.(x), border12(ms)), lc=:black, aspectratio=1,
	title="the phase space",
	xlab=L"m(\Lambda_c^+ K^+)\,\,(\mathrm{MeV})",
	ylab=L"m(\overline{\Lambda}_c^- K^+)\,\,(\mathrm{MeV})")

# ╔═╡ bebec78c-ccd6-4336-9a05-67bb8f2a985c
md"""
The limits of the momentum of $K$ in the B rest frame are determined
by the mass of $\Lambda_c^+ \overline{\Lambda}_c^-$ system since it is two object system, $B\to (\Lambda_c^+ \overline{\Lambda}_c^-) K^+$
"""

# ╔═╡ 9abc515b-d066-478d-8f1b-32bffadab34c
pKofσ(mΛcΛc²) = sqrt(KallenFact(mB^2, mK^2, mΛcΛc²)) /(2mB)

# ╔═╡ 28244da0-f195-4e5e-afb4-2f5c42e52f5d
pKmax = pKofσ.(lims3(ms))[1]

# ╔═╡ fd14cabb-2e77-4a7d-be35-f6c288e51160
md"""
##### The limits are 0 < $P_K$ < $(round.(1000 .* pKmax)) MeV
"""

# ╔═╡ 445d81f9-d336-40a6-9df2-425b6806255d
phsp = flatDalitzPlotSample(ms; Nev=1_000_000) ;

# ╔═╡ 0b66d08d-520b-4110-8850-e6091f04e3a1
let
	xv = pKofσ.(getproperty.(phsp, :σ3))
	stephist(xv,
		xlab=L"P_{B\,\,\mathrm{rest}}(K)\,\,(\mathrm{MeV})",
		title=latexstring("phase space distribution of \$P_{B\\,\\,\\mathrm{ rest}}(K)\$"))
end

# ╔═╡ 4d7a45c0-4b77-4c48-8815-5aeca22bf4e5
md"""
### Check using energy relation
"""

# ╔═╡ 0464c6e1-e2e7-4268-a532-5d23d5783e82
eK(mΛcΛc²) = (mB^2 + mK^2 - mΛcΛc²) /(2mB) # mΛcΛc²(eK) = mB^2 + mK^2 - 2eK*mB

# ╔═╡ c4914e38-4bd4-4b3d-b3f3-860e547a6e1a
pKofe(eK) = sqrt(eK^2-mK^2)

# ╔═╡ 5eeb2014-1e63-4956-a60d-b21a076560c2
pKofe.(eK.(lims3(ms))) # same numbers as above

# ╔═╡ 03bf0964-e32a-48d5-8f16-43bb0892762a
md"""
## Inefficiency
"""

# ╔═╡ 249d37cd-6abf-4f84-91c2-307b4f7441a8
ϵP(p; σ=50e-3) = 1/(1+exp(-p/σ)) ;

# ╔═╡ b42b89f8-7c7a-4de0-9914-a6dd95d2735f
begin
	plot(ϵP, 0, pKmax,
		ylab=L"\varepsilon",
		xlab=L"P_{B\,\,\mathrm{rest}}(K)\,\,(\mathrm{MeV})")
	savefig("assumed_e_pK.pdf")
	plot!()
end

# ╔═╡ 4f0e98ab-0070-49ee-97d0-3c1a293922ab
ϵ(σs) = ϵP(pKofσ(σs[3])) ;

# ╔═╡ b2066aab-7281-4a61-9fd1-b7862eea9fa4
let
	plot(ϵ, ms, clim=(0,1), aspectratio=1,
		title="\$p_K\$ related (in)efficiency",
		xlab=L"m(\Lambda_c^+ K^+)\,\,(\mathrm{MeV})",
		ylab=L"m(\overline{\Lambda}_c^- K^+)\,\,(\mathrm{MeV})")
	savefig("e_on_dalitz.pdf")
	plot!()
end

# ╔═╡ 61132c50-1336-4694-8182-e6ca68e8ade6
let
	plot(e1->project1(ϵ, e1^2, ms) / project1(σs->1.0, e1^2, ms),
		sqrt.(lims1(ms))..., ylim=(0,1),
		xlab=L"m(\overline{\Lambda}_c^- K^+)\,\,(\mathrm{MeV})",
		ylab=L"\epsilon_\mathrm{ph.sp.}")
	savefig("projected1_e.pdf")
	plot!()
end

# ╔═╡ Cell order:
# ╟─9e8a48f6-9e5b-4907-81cd-cabc1ac28213
# ╠═1125f9e0-d4ea-11ec-169b-d754f5700778
# ╠═8b888bf9-7a8f-45e8-b268-41e9fdd8901e
# ╠═f9e829e0-6b08-403d-ba18-697ac51edfad
# ╠═adfc0267-f6e6-412d-80db-3ee0c73bc32b
# ╟─0f28d98f-b7ec-4f98-abd7-017b5b5e4641
# ╟─e2f48768-b995-4260-9e19-a1764af3082e
# ╟─a802e131-2c78-4c3f-9d9f-49aa45f261ac
# ╟─bebec78c-ccd6-4336-9a05-67bb8f2a985c
# ╠═9abc515b-d066-478d-8f1b-32bffadab34c
# ╟─28244da0-f195-4e5e-afb4-2f5c42e52f5d
# ╟─fd14cabb-2e77-4a7d-be35-f6c288e51160
# ╠═445d81f9-d336-40a6-9df2-425b6806255d
# ╟─0b66d08d-520b-4110-8850-e6091f04e3a1
# ╟─4d7a45c0-4b77-4c48-8815-5aeca22bf4e5
# ╠═0464c6e1-e2e7-4268-a532-5d23d5783e82
# ╠═c4914e38-4bd4-4b3d-b3f3-860e547a6e1a
# ╠═5eeb2014-1e63-4956-a60d-b21a076560c2
# ╟─03bf0964-e32a-48d5-8f16-43bb0892762a
# ╠═249d37cd-6abf-4f84-91c2-307b4f7441a8
# ╠═b42b89f8-7c7a-4de0-9914-a6dd95d2735f
# ╠═4f0e98ab-0070-49ee-97d0-3c1a293922ab
# ╟─b2066aab-7281-4a61-9fd1-b7862eea9fa4
# ╟─61132c50-1336-4694-8182-e6ca68e8ade6
