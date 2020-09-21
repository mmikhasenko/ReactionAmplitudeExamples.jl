using ThreeBodyDecay
using PartialWaveFunctions
using Parameters
# 
# 
ms = ThreeBodyMasses(1,1,1; m0=5.0)
σs = randomPoint(ms)
# 
α1,cosβ1x,γ1 = (rand(3) .* 2 .- 1) .* [π,1,π]
σ1,cosθ23x = σs.σ1,cosθ23(σs,ms^2)
# 
τ1 = (σ1, cosβ1x, α1, cosθ23x, γ1)
#
halfintegerspin = (two_j0=1, two_λ0=1, two_s=1, two_λ=1, two_τ=1)
integerspin     = (two_j0=2, two_λ0=2, two_s=2, two_λ=2, two_τ=2)

function check_chain2_1(τ1, ms, pars)
    @unpack two_j0, two_λ0, two_s, two_λ, two_τ = pars
    # 
    τ3 = change_basis_3from1(τ1, ms)
    τ2 = change_basis_2from3(τ3, ms)
    σs = Invariants(σ1=τ1[1], σ2=τ2[1], σ3=τ3[1])
    #
    (σ2, cosβ2x, α2, cosθ31x, γ2) = τ2
    # chain 2
    DD = wignerD_doublearg(two_j0,two_λ0,two_τ,α2,cosβ2x,γ2)
    # 
    Dwd = sum(wignerD_doublearg(two_j0,two_λ0,two_ν,α1,cosβ1x,γ1) * 
              phase(two_ν-two_τ) *
              wignerd_doublearg(two_j0,two_ν,two_τ,cosθhat12(σs,ms^2)) for two_ν in -two_j0:2:two_j0)
    # 
    DD ≈ Dwd, DD, Dwd
end

function check_chain3_1(τ1, ms, pars)
    @unpack two_j0, two_λ0, two_s, two_λ, two_τ = pars
    # 
    τ3 = change_basis_3from1(τ1, ms)
    σs = Invariants(ms; σ1=τ1[1], σ3=τ3[1])
    #
    (σ3, cosβ3x, α3, cosθ12x, γ3) = τ3
    # chain 3
    DD = wignerD_doublearg(two_j0,two_λ0,two_τ,α3,cosβ3x,γ3)
    # 
    Dwd = sum(wignerD_doublearg(two_j0,two_λ0,two_ν,α1,cosβ1x,γ1) * 
              wignerd_doublearg(two_j0,two_ν,two_τ,cosθhat31(σs,ms^2)) for two_ν in -two_j0:2:two_j0)
    # 
    DD ≈ Dwd, DD, Dwd
end
# 
check_chain2_1(τ1, ms, integerspin) # true
check_chain2_1(τ1, ms, halfintegerspin) # false
#
check_chain3_1(τ1, ms, integerspin) # true
check_chain3_1(τ1, ms, halfintegerspin) # true
# 
