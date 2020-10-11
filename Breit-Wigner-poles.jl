### A Pluto.jl notebook ###
# v0.11.2

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

# ╔═╡ 3b6671a0-d32d-11ea-2ea7-25f28d17f69b
md"
## Breit-Wigner amplitude with energy-dependent width
The amplitude is defined as follows:

``
A(s) = \frac{1}{m_R^2-s-im_R\Gamma(s)},
``

where the width incorporates the energy dependence from the phase space available for the decay:

``
\Gamma(s) = \frac{g^2}{16\pi s}\sqrt{(s-(m_1+m_2)^2)(s-(m_1+m_2)^2)},
``
"

# ╔═╡ 91561120-d32c-11ea-1551-a3d6958590ae
using Plots

# ╔═╡ a6a42b20-d32c-11ea-172f-79bd3b94cb58
theme(:wong)

# ╔═╡ a2f9ef30-d329-11ea-28ea-65e3bb0c5e8f
function BW(s; mR, m1, m2, g, ϕ = 0.0)
	1.0/(mR^2-s-1.0im*g^2*
		sqrt(cis(ϕ)*(s-(m1+m2)^2))*sqrt(s-(m1-m2)^2)*cis(-ϕ/2.0) / (16π*s))
end

# ╔═╡ 10815fb2-d32b-11ea-2791-0de8578c4516
@bind g html"<input type=range min=1.0 max=6.0 step=0.1>"

# ╔═╡ bd9b8fae-d329-11ea-1049-392425fc2bf8
@bind mR html"<input type=range min=0.5 max=1.5 step=0.1>"

# ╔═╡ ed9ada90-d329-11ea-3db4-a9516b50f33a
begin
	const mρ = 0.77
	const mπ = 0.139
	const ma₁ = 1.26
end;

# ╔═╡ b7fc3320-d329-11ea-1afd-7db1872cdcea
begin
    amp(s) = BW(s; mR = ma₁, m1 = mρ, m2 =mπ, g = g)
    plot(x->abs2(amp(x+1e-3im)), 0.91^2, 3.0, lab="a₁→ρπ")
    amp(s) = BW(s; mR = mR, m1 = mρ, m2 =mπ, g = g)
    plot!(x->abs2(amp(x+1e-3im)), 0.91^2, 3.0, lab="R→ρπ")
end

# ╔═╡ e8f4d21e-d329-11ea-3545-65d07d5b73e1
let sxv = range(-0.1, 3.0, length=100), syv = range(-1.0, 0.3, length=100)
    cal = [amp(sx+1.0im*sy) for sy in syv, sx in sxv]
    contour(sxv, syv, log.(abs.(cal)), levels=30, xlab="Re s", ylab="Im s")
	hline!([0.0], l=(2,:red), lab="", colorbar=false,
		ann=(2.5, 0, text("Real axis", :red, 11, :bottom)))
end

# ╔═╡ Cell order:
# ╟─3b6671a0-d32d-11ea-2ea7-25f28d17f69b
# ╠═91561120-d32c-11ea-1551-a3d6958590ae
# ╠═a6a42b20-d32c-11ea-172f-79bd3b94cb58
# ╠═a2f9ef30-d329-11ea-28ea-65e3bb0c5e8f
# ╠═10815fb2-d32b-11ea-2791-0de8578c4516
# ╠═bd9b8fae-d329-11ea-1049-392425fc2bf8
# ╠═ed9ada90-d329-11ea-3db4-a9516b50f33a
# ╠═b7fc3320-d329-11ea-1afd-7db1872cdcea
# ╠═e8f4d21e-d329-11ea-3545-65d07d5b73e1
