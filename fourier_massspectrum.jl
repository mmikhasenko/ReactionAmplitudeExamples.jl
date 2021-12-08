### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 55c22390-a272-11eb-36d9-5bff7cbf9da9
using FFTW

# ╔═╡ 9500f160-c0bb-49be-9d69-c5ea77e92c03
using Plots

# ╔═╡ 6ee281b5-828a-436d-8831-8a0b671a00ad
amplitude(s; m=5, Γ=0.05) = 1/(m^2-s-1im*m*Γ)
# amplitude(s; m=5, Γ=0.5) = 1/(m-sqrt(s)-1im*Γ/2)

# ╔═╡ 16146b32-c8b4-44b3-a72c-e5e4bdffdbee
ev = range(1,11,length=200)

# ╔═╡ 50adea0e-853c-4d9d-a91d-3c76f8bfe58c
plot(e->abs2(amplitude(e^2)), ev)

# ╔═╡ 027a69a5-89c4-4976-b217-6f9803afd8da
ampv = map(e->amplitude(e^2), ev)

# ╔═╡ d1e2cec3-c817-4d84-930d-e90ad135c5f5
image = fft(ampv)[1:100]

# ╔═╡ c397126a-ebed-47b7-be51-3e9729bd997b
plot(log.(abs2.(image)))

# ╔═╡ Cell order:
# ╠═55c22390-a272-11eb-36d9-5bff7cbf9da9
# ╠═9500f160-c0bb-49be-9d69-c5ea77e92c03
# ╠═6ee281b5-828a-436d-8831-8a0b671a00ad
# ╠═16146b32-c8b4-44b3-a72c-e5e4bdffdbee
# ╠═d1e2cec3-c817-4d84-930d-e90ad135c5f5
# ╠═c397126a-ebed-47b7-be51-3e9729bd997b
# ╠═50adea0e-853c-4d9d-a91d-3c76f8bfe58c
# ╠═027a69a5-89c4-4976-b217-6f9803afd8da
