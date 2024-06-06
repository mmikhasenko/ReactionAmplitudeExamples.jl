using ThreeBodyDecay

Lb2JpK = ThreeBodySystem(5.62,3.09,0.938,0.49367)

using Plots
using Interact
using LaTeXStrings

function plot_dalitz_with_projections(f=(σ3,σ1)->1.0)
    # layout = @layout [a{0.8h}; grid(1,2)]
    layout = @layout [a{0.65w,0.7h} b; c _]
    plot(layout=layout, size=(800,600), link=:both)
    #
    σ1v = LinRange(Lb2JpK.mthsq[1],Lb2JpK.sthsq[1], 202)
    σ3v = LinRange(Lb2JpK.mthsq[3],Lb2JpK.sthsq[3], 200)
    cal = [Kibble31(σ3,σ1,Lb2JpK) < 0 ? f(σ3,σ1) : NaN for σ3 in σ3v, σ1 in σ1v]
    heatmap!(σ1v, σ3v, cal, c=:viridis, colorbar=false, sp=1,
        ylab="m[J/psi p] (GeV^2)", xlab="m[p K] (GeV^2)")
    #
    calz = map(z->isnan(z) ? 0.0 : z, cal)
    plot!(sum(calz, dims=2)[:,1], σ3v, lab="", xaxis=false, l=(2,:black), sp=2)
    plot!(σ1v, sum(calz, dims=1)[1,:], lab="", yaxis=false, l=(2,:black), sp=3)
end

plot_dalitz_with_projections()

function A(σ3,σ1, CS, Cs)
    σ = [σ1, 0.0, σ3]
    sum(c*amp(σ[ch[1]], ch[2])  for (ch,c) in zip(CS,Cs))
end

I(σ3,σ1, CS, Cs) = abs2(A(σ3,σ1, CS, Cs))

@manipulate for mΛ in LinRange(1.4,1.9,30), ΓΛ in LinRange(0.02,0.15,30)
    plot_dalitz_with_projections((σ3,σ1)->I(σ3,σ1,
        [(1,BreitWigner(mΛ,ΓΛ))],
        [1.0]))
end

Λ1405  = BreitWigner(1.405,   0.090)
Λ1520  = BreitWigner(1.5195,  0.0156)
Λ1690  = BreitWigner(1.685,   0.050)
Λ1810  = BreitWigner(1.80,    0.090)

model = [(1,Λ1405),
 (1,Λ1520),
 (1,BreitWigner(1.6,0.2)),
 (1,Λ1690),
 (1,Λ1810),
 (3,BreitWigner(4.45,0.06))]

# many Lambdas
@manipulate for c1 in LinRange(0,2,30),
    c2 in LinRange(0,2,30),
    c3 in LinRange(0,2,30),
    c4 in LinRange(0,2,30),
    c5 in LinRange(0,2,30)
    #
    plot_dalitz_with_projections((σ3,σ1)->I(σ3,σ1,
        model,
        [c1, c2, c3, c4, c5, 0.0im]))
end

# similar to OK
plot_dalitz_with_projections((σ3,σ1)->I(σ3,σ1,
    model,
    [0.9, 1.0, 0.7, 0.2, 0.3, 0.0im]))

# similar to OK
plot_dalitz_with_projections((σ3,σ1)->I(σ3,σ1,
    model,
    [0.9, 1.0, 0.7, 0.2, 0.3, 0.2im]))
