
using PartialWaveFunctions
using ThreeBodyDecay
using Parameters
using Plots

theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right)

const jω = 1
const λω = 1

const mω  = 0.78266
const mπ⁻ = 0.13957039
const mπ⁰ = 0.1349768
# 
const mb1 = 1.2295 # ± 3.2
const Γb1 = 0.142 # ± 0.009
#
# const mX = 1.7 # GeV
const ms = ThreeBodyMasses(mω, mπ⁻, mπ⁰; m0=1.7)

#############################################################

struct minusone end
import Base:^
^(x::minusone, n::Int) = isodd(n) ? -1 : 1
macro t_str(s::String)
    minusone()
end

function ϵWignerD(j0,ϵP,M,λ, ϕ1,cosθ2,ϕ23)
    M < 0 && return 0.0im
    W⁺ = wignerD(j0, M, λ, ϕ1, cosθ2, ϕ23)
    W⁻ = wignerD(j0,-M, λ, ϕ1, cosθ2, ϕ23)
    # 
    return W⁺ + ϵP * t"-1"^(j0-M) * W⁻
end

# ϵWignerD(1,1,1,1, 0.3,0.3,0.3)


H(j1,λ1,j2,λ2,j,l,s) = 
    clebschgordan(j1,λ1,j2,λ2,s,λ1+λ2) *
        clebschgordan(l,0,s,λ1+λ2,j,λ1+λ2)
# HH(j1,j2,j3,j,l,L,j0) =
#     clebschgordan(j1,λ1,j2,λ2,s,λ1+λ2) *
#         clebschgordan(l,0,s,λ1+λ2,j,λ1+λ2)
#     clebschgordan(j,τ,j3,λ3,S,τ+λ3) *
#         clebschgordan(L,0,S,τ+λ3,j0,τ+λ3)


#############################################################

# with respect to the chain 1
cosζ_for0_sign(k, σs, msq) = cosθhatk1(k, σs, msq), (k==2 ? :ACW : :CW)

cosζ_fori_sign(k, σs, msq) = cosζ(k,σs, msq), (k==2 ? :ACW : :CW)
cosζ_forj_sign(k, σs, msq) = cosζ(k,σs, msq), (k==2 ? :ACW : :CW)
cosζ_fork_sign(k, σs, msq) = cosζ(k,σs, msq), (k==2 ? :ACW : :CW)

abstract type IndexWignerRotation end
struct TrivialWR <: IndexWignerRotation end
struct HatWR{S} <: IndexWignerRotation
    k::Int
end
struct ZetaRepWR{T,S} <: IndexWignerRotation
    k::Int
end
struct ZetaAllWR{S} <: IndexWignerRotation
    k::Int
end

ispositive(wr::Trivial) = true
ispositive(wr::HatWR{S}) where S = S
ispositive(wr::ZetaRepWR{T,S}) where {T,S} = S
ispositive(wr::ZetaAllWR{S}) where S = S

# ind("^0_{1(1)}") -> Trivial()
# # 
# ind("_{1(2)}^0") -> Hat{true}(1)
# ind("_{2(1)}^0") -> Hat{false}(1)
# # 
# ind("_{2(1)}^1") -> ZetaRep{:S,true}(1)
# ind("_{2(1)}^2") -> ZetaRep{:D,true}(2)
# #
# ind("_{1(2)}^3") -> ZetaAll{true}(3)
# ind("_{2(1)}^3") -> ZetaAll{false}(3)

seq(i,j) = (j-i) ∈ (1,-2)
function wr(system_a, reference_b, particle_c=0)
    system_a == reference_b && return TrivialWR()
    S = seq(system_a, reference_b)
    A,B = S ? (system_a, reference_b) : (reference_b, system_a)
    # 
    particle_c == 0 && return Hat{S}(B)
    #
    particle_c ∉ (system_a, reference_b) && return ZetaAll{S}(particle_c)
    #
    T = (particle_c == reference_b) ? :S : :D
    return ZetaRep{T,!(S)}(particle_c)
end

wr(3,1)
wr(1,2)
# 
wr(1,1,2)
wr(1,2,2)
wr(1,2,1)
wr(1,3,2)
wr(2,1,2)
wr(2,1,1)
wr(1,2,3)

struct Chain{T}
    k::Int
    lineshape::T
    jR::Int
    # 
    j0::Int
    P::Int
    M::Int
    ϵ::Int
    # 
    L::Int
    S::Int
    l::Int
    s::Int
end

wignerd_sign(j,λ1,λ2, cosθ, ispositive) =
    (ispositive ? 1 : x"-1"^*(λ1-λ2)) *
        wignerd(j,λ1,λ2, cosθ)

"""
Decay matrix element for  X → ω π⁻ π⁰
"""
function O(chain::Chain, λ, τ, v)
    @unpack σ1, σ2 = v
    @unpack lineshape, k = chain
    @unpack j0, jR = chain
    i,j = ij_from_k(k)
    #
    σs = Invariants(ms; σ1, σ2)
    js = (jω,0,0)
    λs = ( τ,0,0)
    #
    cosθ = cosθij(k, σs, ms^2)
    #
    w0 = wr(k,1,0); cosζ0 = cosζ(w0, σs, ms^2)
    wi = wr(k,1,i); cosζi = cosζ(wi, σs, ms^2)
    wj = wr(k,1,j); cosζj = cosζ(wj, σs, ms^2)
    wk = wr(k,1,k); cosζk = cosζ(wk, σs, ms^2)
    #
    angular = sum(
        wignerd_sign(j0,λ,ν-λs′[k], cosζ0, ispositive(w0)) *
        #
            H(j,τ,js[k],λs′[k],jR,L,S) * x"-1"^(js[k]-λs′[k]) *
        wignerd(j, ν, τ, cosθ) *
            H(js[i],λs′[i],js[j],λs′[j],jR,l,s) * x"-1"^(js[j]-λs′[j]) *
        #
        wignerd_sign(js[i],λs′[i],λs[i], cosζi, ispositive(wi)) *
        wignerd_sign(js[j],λs′[j],λs[j], cosζj, ispositive(wj)) *
        wignerd_sign(js[k],λs′[k],λs[k], cosζk, ispositive(wk))
# 
                for ν in -jR:jR,
                    λs′ in Iterators.product(-js[1]:js[1],
                                             -js[2]:js[2],
                                             -js[3]:js[3]))
    # 
    bw = lineshape(σs[k])
    return angular*bw
end

function amplitude(model::AbstractVector{Chain{T}} where T, v)
    @unpack ϕ_GJ, cosθ_GJ, ϕ_H = v
    @unpack ϕ_ρ, cosθ_ρ, ϕ_π = v
    # 
    a = 0.0im
    for chain in model
        @unpack j0, ϵ, P, M = chain
        a += sum(
            ϵWignerD(j0,ϵ*P,M,λ, ϕ_GJ, cosθ_GJ, ϕ_H)' *
                O(chain,λ,τ,v) *
                wignerD(jω,τ,λω,ϕ_ρ,cosθ_ρ,ϕ_π)'# *
                # Fω(σ1ω,σ2ω)
                    for λ in -j0:j0, τ in -jω:jω)
    end
    return a
end


const ρA = let (mρ,Γρ) = (0.77, 0.15)
    σ->BW(σ,mρ,Γρ)
end


c = Chain(1,ρA,1, 1, -1, 1, 1)

model0 = [c]

typeof(model0) <: Vector{Chain{T}} where T

v0 = let 
    ϕ_GJ, cosθ_GJ, ϕ_H = (0.3,0.3,0.3)
    ϕ_ρ, cosθ_ρ, ϕ_π = (0.3,0.3,0.3)
    @unpack σ1, σ2 = randomPoint(ms)
    (; ϕ_GJ, cosθ_GJ, ϕ_H, σ1, σ2, ϕ_ρ, cosθ_ρ, ϕ_π)
end

function vDalitz(σ1, σ2; v0)
    @unpack ϕ_GJ, cosθ_GJ, ϕ_H, ϕ_ρ, cosθ_ρ, ϕ_π = v0
    (; ϕ_GJ, cosθ_GJ, ϕ_H, σ1, σ2, ϕ_ρ, cosθ_ρ, ϕ_π)
end

intensity(σs; model, v0) = abs2(amplitude(model, vDalitz(σs.σ1, σs.σ2; v0)))

plot(ms, σs->intensity(σs; model=model0, v0))


