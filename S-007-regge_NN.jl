using Plots
using Flux
using Flux: @epochs

const slims = (2.0,20.0);
const mμ = 0.134; const mμsq = mμ^2;

t_sz(s,z) = -(1-z)*(s-4mμsq)/2
u_sz(s,z) = -(1+z)*(s-4mμsq)/2
ampl_sz(s,z) = s^(1/2 + t(s,z)) + s^(1/2 + u(s,z))

const amax = ampl_sz(slims[2],1)

# plot the model
p1 = let
  zv = range(-1,1,length=100)
  sv = range(slims...,length=100)
  calv = ampl_sz.(sv',zv) ./ amax
  heatmap(sv,zv,calv, xlab="s (GeV^2)", ylab="cos(theta)", c=:viridis, clim=(0,1))
end

# data set
dt = let N = 100
  zipx = zip((slims[1]) .+ (slims[2]-slims[1]) .* rand(N), 2 .* rand(N) .- 1)
  xs = [[s,t_sz(s,z)] for (s,z) in zipx]
  ys = [[ampl_sz(s,z)/amax] for (s,z) in zipx]
  zip(xs, ys)
end

# model architecture
model = Chain(
  Dense(2, 4, σ),
  Dense(4, 3, σ),
  Dense(3, 1, σ));
#
loss(x, y) = Flux.mse(model(x), y)
#
opt = ADAM(0.001)
ps = Flux.params(model)
# teach 1000 times with the same data set
@epochs 2000 Flux.train!(loss, ps, dt, opt) #, cb = () -> println("training, $()")

#  plot the model
p2 = let
  zv = range(-1,1,length=100)
  sv = range(slims...,length=100)
  calv = [model([s,t_sz(s,z)])[1] for (z,s) in Iterators.product(zv,sv)]
  heatmap(sv,zv,calv, xlab="s (GeV^2)", ylab="cos(theta)", c=:viridis, clim=(0,1))
end

plot(p1, p2, layout=grid(1,2), size=(1000,350))
