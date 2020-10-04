using ThreeBodyDecay
using PartialWaveFunctions
using LinearAlgebra
using Parameters
using Test
using TypedTables
using Plots


# 
# 
const ms = ThreeBodyMasses(1,1,1; m0=5.0)
#
halfintegerspin = (two_j0=1, two_λ0=1, two_s=1, two_λ=1, two_τ=1)
integerspin     = (two_j0=2, two_λ0=2, two_s=2, two_λ=2, two_τ=2)

function check_chain2_1(τ1, pars, ms=ms)
    @unpack two_j0, two_λ0, two_s, two_λ, two_τ = pars
    (σ1, cosβ1x, α1, cosθ23x, γ1) = τ1
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

function check_chain3_1(τ1, pars, ms=ms)
    @unpack two_j0, two_λ0, two_s, two_λ, two_τ = pars
    (σ1, cosβ1x, α1, cosθ23x, γ1) = τ1
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
# dynamic variables
function randτ1_angles(σs, ms=ms)
    σ1,cosθ23x = σs.σ1,cosθ23(σs,ms^2)
    # polarization angles
    α1,cosβ1x,γ1 = (rand(3) .* 2 .- 1) .* [π,1,π]
    # complete set
    τ1 = (σ1, cosβ1x, α1, cosθ23x, γ1)
    return τ1
end

function randτ1(ms=ms)
    σs = randomPoint(ms)
    randτ1_angles(σs, ms)
end

function rand_dalitz(N, ms=ms)
    σsv = flatDalitzPlotSample(ms, Nev = N)
    return randτ1_angles.(σsv)
end

#
check_chain2_1(randτ1(), integerspin) # true
check_chain2_1(randτ1(), halfintegerspin) # sometimes false
#
check_chain3_1(randτ1(), integerspin) # true
check_chain3_1(randτ1(), halfintegerspin) # sometimes false
#

#
Ry(z) = [wignerd_doublearg(1,i,j,z) for i in [1,-1], j in [1,-1]]
Rz(ϕ) = Diagonal(cis.(ϕ/2 .* [-1,1]))

@test Rz(-0.5) ≈ inv(Rz(0.5))

function R1invRk(τ1,τk)
    ϕk,cosθk,ϕij = τk[3],τk[2],τk[5]
    ϕ1,cosθ1,ϕ23 = τ1[3],τ1[2],τ1[5]
    
    Rk = Rz(ϕk)*Ry(cosθk)*Rz(ϕij)
    R1 = Rz(ϕ1)*Ry(cosθ1)*Rz(ϕ23)
    #
    return inv(R1)*Rk        
end
function R1invR3(τ1, ms=ms)
    τ3 = change_basis_3from1(τ1, ms)
    return R1invRk(τ1,τ3)
end

# [cosθ/2 -sinθ/2
#  sinθ/2  cosθ/2]
function mismatch(τ1, ms=ms)
    M = real(R1invR3(τ1, ms))
    θ′ = 2*atan(M[2,1], M[1,1]) / π * 180
end
# cosθhat31
function mismatch_invariants(τ1, ms=ms)
    τ3 = change_basis_3from1(τ1, ms)
    σs = Invariants(ms; σ1=τ1[1], σ3=τ3[1])
    _cosθhat31 = cosθhat31(σs,ms^2)
    θhat31 = acos(_cosθhat31) / π * 180
    return θhat31
end

# MC
ll = let
    τ1v = rand_dalitz(10_000)
    τ3v = change_basis_3from1.(τ1v, Ref(ms))
    θ′v = mismatch.(τ1v)
    θv  = mismatch_invariants.(τ1v)
    fv = abs.(θv .- θ′v) .< 1e-5

    [(NamedTuple{(:σ1, :cosβ1, :α1, :cosθ1, :γ1)}(τ1)...,
      NamedTuple{(:σ3, :cosβ3, :α3, :cosθ3, :γ3)}(τ3)...,
      θw_so2 = θ′, θw_inv = θ, zero=v) for (τ1,τ3,v,θ,θ′) in zip(τ1v,τ3v,fv,θv,θ′v)]
end
tll = Table(ll)

let
    scatter(tll.α1[iszero.(tll.zero)], tll.α3[iszero.(tll.zero)], lab="need +2π correction")
    scatter!(tll.α1[tll.zero], tll.α3[tll.zero], lab="does no need correction")
    plot!(xlab="ϕ₁", ylab="ϕ₃", size=(450,400))
end

nocorrection(x) = abs(x.α1-x.α3) < π
@test filter(nocorrection, tll) == filter(x->x.zero, tll)

let
    ttl_s = filter(x->abs(x.α1+0.3*x.α3) < 0.1π, tll)
    # scatter!(ttl_s.α1, ttl_s.α3, lab="")
    # 
    scatter(ttl_s.α3, ttl_s.θw_inv, lab="")
    scatter!(ttl_s.α3, ttl_s.θw_so2, lab="")
end
