using AlgebraPDF
using Plots

const mπ = 0.15
const mρ = 0.77
const Γρ = 0.15
# 
k(s) = sqrt(s/4-mπ^2)
ρ(s) = 2k(s)/(8π*√s)
# 
const gsq = 2mρ*Γρ/ρ(mρ^2)
const a = gsq/(8π*mρ*(mρ^2-4mπ^2)) # 1/GeV
const r = -4mρ/gsq*8π*2+.5 # 1/GeV
#
const UnitGeVxfm = 0.2
# 
a*UnitGeVxfm # 0.13 fm
r*UnitGeVxfm # -4.9 fm

1-ca/(2r)


fBW  = Normalized(FunctionWithParameters((s;p)->1/abs2(mρ^2-s-1im*mρ*Γρ), ∅), (0.1, 1.1))
fBWg = Normalized(FunctionWithParameters((s;p)->1/abs2(mρ^2-s-1im*gsq*ρ(s)/2), ∅), (0.1, 1.1))
fSL  = Normalized(FunctionWithParameters((s;p)->1/abs2(1/a+r/2*k(s)^2-1im*k(s)), ∅), (0.1, 1.1))
#

plot()
plot!(fBW)
plot!(fBWg)
plot!(fSL)

Z(a,r) = 1-sqrt(1/(1+2*abs(r/a)))

Z(a,r)