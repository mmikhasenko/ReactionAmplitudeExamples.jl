using QuadGK
using Plots
theme(:wong)
using Test

"""
    Isc(m1sq, m2sq, m3sq, M1sq, M2sq, M3sq)

    Scalar triangle loop function
              o --- M₂
          m₁ / \\
            /   \\ m₃
    M₃ --- o --- o --- M₁
              m₂
"""
function Isc(m1sq, m2sq, m3sq, M1sq, M2sq, M3sq)
    function integrand(x)
        A = M1sq;
        B = m2sq + x*M2sq - (1-x)*M1sq - x*M3sq-m3sq;
        C = x*m1sq + (1-x)*m3sq - x*(1-x)*M2sq;
        D = B^2-4*A*C;
        y1, y2 = (-B .- sqrt(D) .* [1.0, -1.0]) ./ (2*A);
        return (log((1-x-y2)/(-y2))-log((1-x-y1)/(-y1)))/(A*(y2-y1))
    end
    integr = quadgk(integrand, 0, 1)[1];
    return integr/(16*π^2);
end


# masses in GeV
const mD0 = 1.86483; const mD0sq = mD0^2
const mDstar0 = 2.00685; const mDstar0sq = mDstar0^2
const ΓDstar = 65e-6;
#
const mπ0 = 0.13498; const mπ0sq = mπ0^2
const mπ = 0.13957; const mπsq = mπ^2

F(s,σ) = Isc(mDstar0sq-1im*mDstar0*ΓDstar, mD0sq, mπsq, σ, mD0sq, s)

# energy in MeV counted from threshold
e2s(e) = (mD0+mDstar0+e*1e-3)^2
e2σ(e) = (mD0+mπ+e*1e-3)^2
# 
@test e2σ(0)<mDstar0sq<e2σ(3)
@test e2s(-4)<(2mD0+mπ)^2<e2s(0)

# check if TS is well seen
let
    ev = range(-1,3,length=100)
    calv = F.(e2s.(ev), e2σ(1))
    plot(ev, [real.(calv) imag.(calv)], lab=["re" "im"])
end

λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
ρ(σ) = sqrt(λ(σ,mD0sq,mπsq))/mDstar0sq/(8π)
# 
const ρDstar0 = ρ(mDstar0sq)
const gsq = 2*mDstar0*ΓDstar / ρDstar0
t(σ) = gsq/(mDstar0sq-σ-1im*gsq*ρ(σ)/2)

a(s,σ) = t(σ)*(1+gsq*F(s,σ))

let s = e2s(1.1)
    plot(xlab="m(D⁰π)", title="corrected lineshape")
    ev = range(0,3,length=300)
    calv = a.(s, e2σ.(ev))
    plot!(ev, [real.(calv) imag.(calv)], lab=["re[a]" "im[a]"], seriescolor=[1 2])
    calv = t.(e2σ.(ev))
    plot!(ev, [real.(calv) imag.(calv)], lab=["re[t]" "im[t]"], seriescolor=[1 2], ls=:dash)
end

# quasi two-body phase space
integrand_t(s,σ) = sqrt(λ(s,σ,mD0sq)*λ(σ,mD0sq,mπsq))/(σ*s)*abs2(t(σ))
integrand_a(s,σ) = sqrt(λ(s,σ,mD0sq)*λ(σ,mD0sq,mπsq))/(σ*s)*abs2(a(s,σ))
ρ3_a(s) = quadgk(σ->integrand_a(s,σ), (mπ+mD0)^2, (sqrt(s)-mD0)^2)[1]
ρ3_t(s) = quadgk(σ->integrand_t(s,σ), (mπ+mD0)^2, (sqrt(s)-mD0)^2)[1]


# effect of the exchange: t vs a
let
    plot(xlab="m(D⁰D⁰π)")
    ev = range(-1,1,length=50).*(5e3*ΓDstar)
    calv_t = ρ3_t.(e2s.(ev))
    plot!(ev, calv_t, lab="ρ₃ bubble")
    # 
    calv_a = ρ3_a.(e2s.(ev))
    plot!(ev, calv_a, lab="ρ₃ bubble+exchange")
    # 
    plot!(ev, calv_a-calv_t, lab="ρ₃ exchange", seriescolor=6)
end
