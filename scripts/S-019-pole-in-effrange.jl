using Plots
using AlgebraPDF
using QuadGK
using Interpolations


theme(:wong2, grid=false, frame=:box)

const mD = 2.0
k(s) = sqrt(s/4-mD^2)

A(s;a,r) = (_k = k(s); 1/(1/a+r/2*_k^2 - 1im*_k))
A_II(s;a,r) = (_k = k(s); 1/(1/a+r/2*_k^2 + 1im*_k))

e2m(e) = 2mD+e

X = FunctionWithParameters(
    (e;p)->A(e2m(e+1e-6im)^2;a=p.a,r=p.r), (a=0.1,r=0.0))
# 
X_II = FunctionWithParameters(
    (e;p)->A_II(e2m(e+1e-6im)^2;a=p.a,r=p.r), (a=0.1,r=0.0))

plot(abs2(X), -0.1, 0.1)
plot(abs2(X_II), -0.1, 0.1)

let a = -1.9, r=0.96
    _X = updatepars(X,(;a,r))
    _X_II = updatepars(X_II,(;a,r))
    plot(size=(300,500), layout=grid(2,1),
        contour(-1:0.01:1, -1:0.01:1, (x,y)->log(abs2(_X(x+1im*y))), colorbar=false),
        contour(-1:0.01:1, -1:0.01:1, (x,y)->log(abs2(_X_II(x+1im*y))), colorbar=false))
end

const m1 = 0.2
const m2 = 0.1
const ΓD = 0.4

quadgk(x->x^2, 1.0, 1.0+1im, 3.0+3im)[1]

λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
k3b(s::Complex) = quadgk(
        # σ->sqrt(λ(s,σ,mD^2))/
        # σ->sqrt((s-(sqrt(σ)+mD)^2))*sqrt((s-(sqrt(σ)-mD)^2))/
        # σ->sqrt(σ-(sqrt(s)+mD)^2)*sqrt(σ-(sqrt(s)-mD)^2)/
        σ->sqrt((sqrt(s)+mD)^2-σ)*sqrt((sqrt(s)-mD)^2-σ)/
            (mD^2-σ-1im*mD*ΓD)/(mD^2-σ+1im*mD*ΓD),
        (m1+m2)^2,
        real((sqrt(s)-mD)^2)+sign(imag(s))*1e-6im,
        (sqrt(s)-mD)^2)[1] * ΓD / 3 /
    sqrt(s)
#
plot( e->real(k3b(e2m(e+1e-6im)^2)), -1, 1)
plot!(e->real(k(e2m(e+1e-6im)^2)), -1, 1)
# 
contour(-2:0.01:1, -1:0.01:1,
    (x,y)->log(abs2(k3b(e2m(x+1im*y)^2))), colorbar=false)
#


const ecut = 2.0
const scut = e2m(ecut)^2
#
const xv = range((m1+m2+mD)^2, scut, length=100)
const yv = real.(k3b.(xv .+ 1e-7im))
const itr = interpolate((xv,), yv, Gridded(Linear()))
# 
const cutratio = yv[end] / k(scut)

k3b(s::Real) = s<scut ? itr(s) : cutratio*k(s)

let
    plot( e->real(k3b(e2m(e+1e-6im)^2)), -1, 4)
    plot!(e->real(k(e2m(e+1e-6im)^2)), -1, 4)
    plot!(e->real(k3b(e2m(e)^2)), -1, 4)
    plot!(e->imag(ik3b_tilde(e2m(e+1e-6im)^2)), -1, 4)
end
# 
ik3b_tilde(s::Complex) = s*sqrt(s)*quadgk(s′->k3b(s′)/sqrt(s′)/s′/(s′-s),(m1+m2+mD)^2, Inf)[1]/π
const ksub = real(ik3b_tilde(e2m(0.0)^2+1e-6im))
ik3b_tilde_sub(s::Complex) = ik3b_tilde(s)-ksub

ik3b_tilde_sub_II(s::Complex) = ik3b_tilde_sub(s)+2im*k3b(s)

@time ik3b_tilde_sub(e2m(0.0)^2+1e-6im)

ik3b_tilde_sub(e2m(1.0)^2-1e-6im)
ik3b_tilde_sub(e2m(1.0)^2+1e-6im)
ik3b_tilde_sub_II(e2m(1.0)^2-1e-6im)

k3b(e2m(1.0)^2)

B(s;a,r) = (_ik = ik3b_tilde_sub(s); 1/(1/a+r/2*k(s)^2 - _ik))
B_II(s;a,r) = (_ik = ik3b_tilde_sub_II(s); 1/(1/a+r/2*k(s   )^2 - _ik))
B_I_II(s;a,r) = imag(s) > 0.0 ? B(s;a,r) : B_II(s;a,r)

let a = -10.0, r=0.55
    contour(size=(400,300),
        -2:0.03:1, -1:0.03:0.3, (x,y)->log(abs2(B_I_II(e2m(x+1im*y)^2;a,r))),
        colorbar=false, levels=50)
    # plot!([-1.7,1],[0,0], l=(2,:red), lab="")
    scatter!([-1.7],[0], m=(4,:red), lab="")
end