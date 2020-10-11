### A Pluto.jl notebook ###
# v0.11.1

using Markdown
using InteractiveUtils

# ╔═╡ 96942290-d326-11ea-38b4-418a0b8da41a
using Plots

# ╔═╡ 9b4bde40-d326-11ea-3bd5-378203569b99
using Flux

# ╔═╡ a0e7432e-d326-11ea-1615-756b95e4a849
using Flux: @epochs

# ╔═╡ aad9cd40-d326-11ea-1d83-eb16bf33b9d5
begin
	const slims = (2.0,20.0);
	const mμ = 0.134;
	const mμsq = mμ^2;
end;

# ╔═╡ b6eb8dd0-d326-11ea-176c-058904eb0a2a
begin
	t_sz(s,z) = -(1-z)*(s-4mμsq)/2
	u_sz(s,z) = -(1+z)*(s-4mμsq)/2
	ampl_sz(s,z) = s^(1/2 + t_sz(s,z)) + s^(1/2 + u_sz(s,z))
end

# ╔═╡ ba7291b0-d326-11ea-2863-d5b44790393c
const amax = ampl_sz(slims[2],1)

# ╔═╡ 845ff2c0-d326-11ea-3b06-9fa9f292ef5e
# plot the model
p1 = let
  zv = range(-1,1,length=100)
  sv = range(slims...,length=100)
  calv = ampl_sz.(sv',zv) ./ amax
  heatmap(sv,zv,calv, xlab="s (GeV^2)", ylab="cos(theta)", c=:viridis, clim=(0,1))
end

# ╔═╡ 4b4d4850-d328-11ea-2839-69f24a7f0418
# data set
dt = let N = 100
  zipx = zip((slims[1]) .+ (slims[2]-slims[1]) .* rand(N), 2 .* rand(N) .- 1)
  xs = [[s,t_sz(s,z)] for (s,z) in zipx]
  ys = [[ampl_sz(s,z)/amax] for (s,z) in zipx]
  zip(xs, ys)
end

# ╔═╡ ff0e3d50-d327-11ea-2b91-9b5afc6a2b0c
# model architecture
model = Chain(
  Dense(2, 4, σ),
  Dense(4, 3, σ),
  Dense(3, 1, σ));

# ╔═╡ 08c6f2b0-d328-11ea-301e-2b030656d147
loss(x, y) = Flux.mse(model(x), y)

# ╔═╡ 0be570c0-d328-11ea-34f4-0729c1a0007a
begin
	opt = ADAM(0.01)
	ps = Flux.params(model)
end;

# ╔═╡ 14a79120-d328-11ea-2988-39eb7d0d047f
@epochs 100 Flux.train!(loss, ps, dt, opt) #, cb = () -> println("training, $()")

# ╔═╡ 196407c0-d328-11ea-3e55-452963d31671
#  plot the model
p2 = let
	zv = range(-1,1,length=100)
	sv = range(slims...,length=100)
 	calv = [model([s,t_sz(s,z)])[1] for (z,s) in Iterators.product(zv,sv)]
 	heatmap(sv,zv,calv, xlab="s (GeV^2)", ylab="cos(theta)", c=:viridis, clim=(0,1))
end

# ╔═╡ Cell order:
# ╠═aad9cd40-d326-11ea-1d83-eb16bf33b9d5
# ╠═b6eb8dd0-d326-11ea-176c-058904eb0a2a
# ╠═ba7291b0-d326-11ea-2863-d5b44790393c
# ╠═845ff2c0-d326-11ea-3b06-9fa9f292ef5e
# ╠═4b4d4850-d328-11ea-2839-69f24a7f0418
# ╠═ff0e3d50-d327-11ea-2b91-9b5afc6a2b0c
# ╠═08c6f2b0-d328-11ea-301e-2b030656d147
# ╠═0be570c0-d328-11ea-34f4-0729c1a0007a
# ╠═14a79120-d328-11ea-2988-39eb7d0d047f
# ╠═196407c0-d328-11ea-3e55-452963d31671
# ╠═96942290-d326-11ea-38b4-418a0b8da41a
# ╠═9b4bde40-d326-11ea-3bd5-378203569b99
# ╠═a0e7432e-d326-11ea-1615-756b95e4a849
