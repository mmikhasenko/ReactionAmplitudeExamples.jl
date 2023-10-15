### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 1e4eda37-22b6-493c-968a-57f66895bf32
using Optim

# ╔═╡ b87c2e95-db4c-4de0-a8a5-8ed545d92d9f
λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x

# ╔═╡ 1c937b1d-ec6b-446a-97a8-d97ab17d8952
begin
	mπ⁺ = 0.139570
	mD⁰, mD⁺ = 1.86483, 1.86965
	mDˣ⁺ = 2.01026
end

# ╔═╡ ee0d7bfc-76cd-41bc-adf9-ae58e6c3f27a
mth = mDˣ⁺+mD⁰

# ╔═╡ 85552381-75d8-462e-9138-db579ccf007b
begin
	δm0 = 281e-6
	m0, Γ0 = mth-δm0, 0.394e-3
end

# ╔═╡ 0149a2a0-a292-11eb-39a8-79ebcc4f4f7f
begin
	p(s) = sqrt(λ(s,mD⁰^2,mπ⁺^2)/(4s))
	ρ(s) = 2*p(s)/sqrt(s)
	Γ(s) = Γ0 * ρ(s)/ρ(m0^2) * p(s)^2/p(m0^2)^2 # P-wave
end;

# ╔═╡ 0af8042e-a887-4d9e-830e-ddc80d8286a8
BW(s) = m0*Γ0/(m0^2-s-1im*m0*Γ(s))

# ╔═╡ c6cd5084-7b28-4164-b49f-44be26ae967c
to_minimize(r,i) = 1/BW((r-1im*i/2)^2) |> abs2

# ╔═╡ fe30c847-bd9e-4107-8472-3c01b9581fd9
md"""
The values of $m_0$ and $Γ_0$ are not exactly the pole position.
Here we check it:
"""

# ╔═╡ e4e1dcc5-ff6e-47a2-b94c-ca8c7b2d8697
to_minimize(m0,Γ0)

# ╔═╡ cc0c9447-96fa-42aa-a365-101b791b9054
fit = optimize(x->to_minimize(x...), [m0,Γ0])

# ╔═╡ 9c371ea8-7a56-4f00-86b9-274dad58b2fb
mΓ = NamedTuple{(:m_pole,:Γ_pole)}(fit.minimizer)

# ╔═╡ 793e7258-7abf-4623-8ecd-e812cba5e2fe
md"""
Here is a check that the denominator vanishes at the pole
"""

# ╔═╡ 7a23c1e9-154c-4b82-bdce-ac83929c2363
to_minimize(mΓ...) < 1e-7

# ╔═╡ b055489f-8b67-4212-bff3-0dfaaf82a1a1
mΓ.m_pole-mth

# ╔═╡ Cell order:
# ╠═1e4eda37-22b6-493c-968a-57f66895bf32
# ╠═b87c2e95-db4c-4de0-a8a5-8ed545d92d9f
# ╠═1c937b1d-ec6b-446a-97a8-d97ab17d8952
# ╠═ee0d7bfc-76cd-41bc-adf9-ae58e6c3f27a
# ╠═85552381-75d8-462e-9138-db579ccf007b
# ╠═0149a2a0-a292-11eb-39a8-79ebcc4f4f7f
# ╠═0af8042e-a887-4d9e-830e-ddc80d8286a8
# ╠═c6cd5084-7b28-4164-b49f-44be26ae967c
# ╟─fe30c847-bd9e-4107-8472-3c01b9581fd9
# ╠═e4e1dcc5-ff6e-47a2-b94c-ca8c7b2d8697
# ╠═cc0c9447-96fa-42aa-a365-101b791b9054
# ╠═9c371ea8-7a56-4f00-86b9-274dad58b2fb
# ╟─793e7258-7abf-4623-8ecd-e812cba5e2fe
# ╠═7a23c1e9-154c-4b82-bdce-ac83929c2363
# ╠═b055489f-8b67-4212-bff3-0dfaaf82a1a1
