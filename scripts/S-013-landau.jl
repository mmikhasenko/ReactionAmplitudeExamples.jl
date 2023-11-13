using QuadGK
using AlgebraPDF
using Plots
using Interpolations


theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right)

integrand(t,x) = exp(-t)*cos(t*x+2*t/π*log(t))
landau(x) = quadgk(t->integrand(t,x), 0, Inf)[1]

@time landau(0)
@time plot(x->landau(x), -2:0.5:15)

const xv = range(-3,30, length=300)
const yv = landau.(xv)



const itr = interpolate((xv,), yv, Gridded(Linear()))
dlandau(x) = itr.knots[1][1]<x<itr.knots[1][end] ? itr(x) : 0.0

plot(dlandau, -10, 40)



import AlgebraPDF:func
struct Landau{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::Landau, x::NumberOrTuple; p=pars(d))
    allp = p+fixedpars(d)
    μ,c = (getproperty(allp, s) for s in keys(d.p))
    dlandau(x-μ/c)
end

d = Landau((μ=0.0,c=1.0))
d(0.0)
@time plot(d, -5, 15)

d1 = Landau(Ext(μ1= 2.0,c1=1.0)) |> d->fixpars(d,(:μ1,:c1))
d2 = Landau(Ext(μ2= 10.5,c2=1.2)) |> d->fixpars(d,(:μ2,:c2))
d3 = Landau(Ext(μ3=-2.5,c3=1.5)) |> d->fixpars(d,(:μ3,:c3))

dsum = d1*(α1=1.1,) + d2*(α2=1.0,) + d3*(α3=1.1,)

plot(dsum, -5, 15)

pars(dsum) #
# (α1 = 1.1, α2 = 1.0, α3 = 1.1, μ1 = 2.0, c1 = 1.0, μ2 = 10.5, c2 = 1.2, μ3 = -2.5, c3 = 1.5)

data = generate(Normalized(dsum, (-5,15)), 10000)



h = Plots.StatsBase.fit(
    Plots.StatsBase.Histogram,
    data, nbins=100)
const xdata = Plots._bin_centers(h.edges[1])
const ydata = h.weights

plot(h, c=1, lab="data")








import AlgebraPDF: func, pars, updatevalueorflag
struct ChiSq{M,T<:NumberOrTuple,V<:Number} <: AbstractFunctionWithParameters
    f::M
    xdata::Vector{T}
    ydata::Vector{V}
end
function func(d::ChiSq, x::NumberOrTuple; p=pars(d))
    sum((ydata - dsum(xdata; p)).^2) #  ./ ydata
end
pars(d::ChiSq, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag(d::ChiSq, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    ChiSq(updatevalueorflag(d.f,s,isfree,v), d.xdata, d.ydata)

# 

χ² = ChiSq(dsum, collect(xdata), ydata)


χ²(0,p2v(χ²))


fit_result = AlgebraPDF.minimize(x->χ²(0,x), p2v(χ²), MigradAndHesse())
fs = AlgebraPDF.fit_summary(fit_result, χ²)

fs.parameters

plot!(updatepars(dsum, fs.parameters), xdata[1], xdata[end], l=(2,:red), lab="fit")
