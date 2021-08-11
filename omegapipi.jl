
using PartialWaveFunctions
using ThreeBodyDecay
using Parameters
using Plots

theme

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

#############################################################

struct Chain{T}
    k::Int
    lineshape::T
    j::Int
    # 
    j0::Int
    P::Int
    M::Int
    ϵ::Int
end

"""
Decay matrix element for  X → ω π⁻ π⁰
"""
function O(chain::Chain, λ, τ, v)
    @unpack σ1, σ2 = v
    @unpack lineshape, k = chain
    @unpack j0, j = chain
    #
    _σs = Invariants(ms; σ1, σ2)
    _cosθhat = cosθhatk1(k,_σs, ms^2)
    _cosθij = cosθij(k, _σs, ms^2)
    _σk = _σs[k]
    # 
    angular = sum(
        wignerd(j0,λ,ν, _cosθhat)*
        wignerd(j, ν,τ, _cosθij ) for ν in -j:j)
    # 
    bw = lineshape(_σk)
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


