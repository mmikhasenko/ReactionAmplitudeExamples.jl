
using PartialWaveFunctions
using ThreeBodyDecay
using DataFrames
using Parameters
using Plots
using AlgebraPDF
using LaTeXStrings
using DelimitedFiles
using LinearAlgebra
import Plots.PlotMeasures.mm


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
const jps3 = (jp(jω,'-'), jp"0-", jp"0-")

#############################################################

ifsmallgiveanan(v; cut=1e-8) = abs(v)< cut ? NaN : v

#############################################################

struct minusone end
import Base:^
^(x::minusone, n::Int) = isodd(n) ? -1 : 1
macro x_str(s::String)
    minusone()
end

function ϵWignerD(j0,ϵP,M,λ, ϕ1,cosθ2,ϕ23)
    M < 0 && return 0.0im
    W⁺ = wignerD(j0, M, λ, ϕ1, cosθ2, ϕ23)
    W⁻ = wignerD(j0,-M, λ, ϕ1, cosθ2, ϕ23)
    # 
    return W⁺ + ϵP * x"-1"^(j0-M) * W⁻
end

wignerd_sign(j,λ1,λ2, cosθ, ispositive) =
    (ispositive ? 1 : x"-1"^*(λ1-λ2)) *
        wignerd(j,λ1,λ2, cosθ)

function pk(k,σs,msq)
    i,j,_ = ijk(k)
    sqrt(ThreeBodyDecay.λ(σs[k],msq[i],msq[j])/(4σs[k]))
end
qk(k,σs,msq) = sqrt(ThreeBodyDecay.λ(msq[4],σs[k],msq[k])/(4msq[4]))


H(j1,λ1,j2,λ2,j,l,s) = 
    clebschgordan(j1,λ1,j2,-λ2,s,λ1-λ2) *
        clebschgordan(l,0,s,λ1-λ2,j,λ1-λ2)
# 

#############################################################

@with_kw struct Chain{T}
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

function nt(c::Chain)
    ks = fieldnames(Chain)
    NamedTuple{ks}(getfield(c,k) for k in ks)
end

import ThreeBodyDecay:jp
jp(c::Chain) = jp(c.j0,c.P==0 ? '-' : '+')


pm(v) = ! iszero(v) ? '+' : '-'
ξn(c::Chain) = c.k==1 ? ["f_0", "\\rho", "f_2", "\\rho_3"][c.jR+1] : ["?", "b_1", "?"][c.jR+1]
πn(c::Chain) = c.k==1 ? "\\omega" : "\\pi"
Lwave(l::Int) = l ≤ 6 ? ['S','P','D','F','G','H','I'][l+1] : string(l)[1]
# 
JPCMϵ(c::Chain) = "$(c.j0)^{$(pm(c.P))+} $(c.M)^{$(pm(c.ϵ))}"
ξL(c::Chain) = "$(ξn(c))$(πn(c))\\,\\,$(Lwave(c.L))"
ξLS(c::Chain) = "\\left[$(ξn(c))$(πn(c))\\right]_{$(c.S)}\\,\\,$(Lwave(c.L))"
# 
JPCMϵξL(c::Chain) = JPCMϵ(c)*" "*ξL(c)

                                                   
#                                  _|            _|  
#  _|_|_|  _|_|      _|_|      _|_|_|    _|_|    _|  
#  _|    _|    _|  _|    _|  _|    _|  _|_|_|_|  _|  
#  _|    _|    _|  _|    _|  _|    _|  _|        _|  
#  _|    _|    _|    _|_|      _|_|_|    _|_|_|  _|  




"""
Decay matrix element for  X → ω π⁻ π⁰
"""
function O(chain::Chain, λ, τ, v)
    @unpack σ1, σ2 = v
    @unpack lineshape, k = chain
    @unpack j0, jR = chain
    @unpack l,s,L,S = chain
    i,j,_ = ijk(k)
    #
    σs = Invariants(ms; σ1, σ2)
    js = getproperty.(jps3,:j)
    λs = ( τ,0,0)
    #
    cosθ = cosθij(k, σs, ms^2)
    p = pk(k,σs,ms^2)
    q = qk(k,σs,ms^2)
    #
    w0 = wr(k,1,0); cosζ0 = cosζ(w0, σs, ms^2)
    wi = wr(k,1,i); cosζi = cosζ(wi, σs, ms^2)
    wj = wr(k,1,j); cosζj = cosζ(wj, σs, ms^2)
    wk = wr(k,1,k); cosζk = cosζ(wk, σs, ms^2)
    #
    angular = sum(
        wignerd_sign(j0,λ,ν-λs′[k], cosζ0, ispositive(w0)) *
        #
        q^L * H(jR,ν,js[k],λs′[k],j0,L,S) * x"-1"^(js[k]-λs′[k]) *
        wignerd(jR, ν, λs′[i]-λs′[j], cosθ) *
        p^L * H(js[i],λs′[i],js[j],λs′[j],jR,l,s) * x"-1"^(js[j]-λs′[j]) *
        #
        wignerd_sign(js[i],λs′[i],λs[i], cosζi, ispositive(wi)) *
        wignerd_sign(js[j],λs′[j],λs[j], cosζj, ispositive(wj)) *
        wignerd_sign(js[k],λs′[k],λs[k], cosζk, ispositive(wk)) *
        1.0
        # 
                for ν in -jR:jR,
                    λs′ in Iterators.product(-js[1]:js[1],
                                             -js[2]:js[2],
                                             -js[3]:js[3]))
    # 
    bw = lineshape(σs[k])
    return angular*bw
end

Osym(chain::Chain, λ, τ, v) = 
    chain.k==1 ? 
        O(chain,λ,τ,v) :
        (O(Chain(chain; k=2),λ,τ,v) + # SIGN TO BE CHECKED!
            O(Chain(chain; k=3),λ,τ,v))

function amplitude10d(model::AbstractVector{Chain{T}} where T, v)
    @unpack ϕ_GJ, cosθ_GJ, ϕ_H = v
    @unpack ϕ_ρ, cosθ_ρ, ϕ_π = v
    # 
    a = 0.0im
    for chain in model
        @unpack j0, ϵ, P, M = chain
        a += sum(
            ϵWignerD(j0,ϵ*P,M,λ, ϕ_GJ, cosθ_GJ, ϕ_H)' *
                Osym(chain, λ, τ, v) *
                wignerD(jω,τ,λω,ϕ_ρ,cosθ_ρ,ϕ_π)'# *
                # Fω(σ1ω,σ2ω)
                    for λ in -j0:j0, τ in -jω:jω)
    end
    return a
end

# test

testchain = Chain(k=1,lineshape=x->BW(x,0.77,0.15),jR=1,
    j0=1,P=1,M=1,ϵ=1,
    L=1,S=2,l=1,s=0)

testmodel = [testchain]

testvars = let 
    ϕ_GJ, cosθ_GJ, ϕ_H = (0.3,0.3,0.3)
    ϕ_ρ, cosθ_ρ, ϕ_π = (0.3,0.3,0.3)
    @unpack σ1, σ2 = randomPoint(ms)
    (; ϕ_GJ, cosθ_GJ, ϕ_H, σ1, σ2, ϕ_ρ, cosθ_ρ, ϕ_π)
end

@assert amplitude10d(testmodel, testvars) != 0.0

function vDalitz(σ1, σ2; v0)
    @unpack ϕ_GJ, cosθ_GJ, ϕ_H, ϕ_ρ, cosθ_ρ, ϕ_π = v0
    (; ϕ_GJ, cosθ_GJ, ϕ_H, σ1, σ2, ϕ_ρ, cosθ_ρ, ϕ_π)
end

intensity10d(σs; model, v0) = abs2(amplitude10d(model, vDalitz(σs.σ1, σs.σ2; v0)))

# plots Dalitz
plot(ms, σs->intensity10d(σs; model=testmodel, v0=testvars))



function intensitydalitz(σs; coherentchains)
    j0 = coherentchains[1].j0
    @unpack σ1, σ2 = σs
    v = (; σ1,σ2)
    # 
    return sum(abs2, 
        sum(Osym(chain,λ,τ,v) for chain in coherentchains)
            for λ in -j0:j0, τ ∈ -jω:jω)
end

@recipe f(chain::Chain) = 
    (ms,σs->intensitydalitz(σs; coherentchains=[chain]))


                                                               
                                                               
#  _|      _|      _|    _|_|_|  _|      _|    _|_|      _|_|_|  
#  _|      _|      _|  _|    _|  _|      _|  _|_|_|_|  _|_|      
#    _|  _|  _|  _|    _|    _|    _|  _|    _|            _|_|  
#      _|      _|        _|_|_|      _|        _|_|_|  _|_|_|    



# Isobars
fρ  = FBreitWigner((mρ=0.7685, Γρ=0.1507))
fρ3 = FBreitWigner((mρ3=1.690, Γρ3=0.190))
ff2 = FBreitWigner((mf2=1.274, Γf2=0.185))
fb1 = FBreitWigner((;mb1, Γb1))

resonances = [
    (n=:ρ,  a = fρ,  jp=jp"1-", k=1),
    (n=:f2, a = ff2, jp=jp"2+", k=1),
    (n=:ρ3, a = fρ3, jp=jp"3-", k=1),
    (n=:b1, a = fb1, jp=jp"1+", k=2),
]

let
    plot()
    for (n,a) in resonances
        plot!(abs2(a) , 0.1, 1.7, lab="$n")
    end
    plot!()
end

function selectwaveset(JP; Mmax=1, Lmax=3, ϵ=1)
    waveset = DataFrame(
        k=Int[], lineshape=AbstractFunctionWithParameters[], jR=Int[],
        j0=Int[], P=Int[], M=Int[], ϵ=Int[],
        # ls=NTuple{2,Int}[], LS=NTuple{2,Int}[])
        l=Int[], s=Int[], L=Int[], S=Int[])
    # 
    for (n,a,jp,k) in resonances
        # 
        for M in 0:JP.j,
            lsLS in possible_lsLS(k, jp, [jps3...,JP])
            
            P = JP.p=='+' ? 1 : -1
            (M==0 && (ϵ*P*x"-1"^JP.j) == -1) && continue
            #
            lsLS.LS[1] > Lmax && continue
            M > Mmax && continue
            # 
            push!(waveset, (; k,lineshape=a,jR=jp.j,
                        j0=JP.j,P,M,ϵ,
                        NamedTuple{(:l,:s)}(lsLS.ls)...,
                        NamedTuple{(:L,:S)}(lsLS.LS)...))
        end
    end
    return waveset
end


function crossedterm(σs, chain1::Chain, chain2::Chain)
    chain1.j0 != chain2.j0 && return 0.0im
    chain1.M  != chain2.M  && return 0.0im
    chain1.P  != chain2.P  && return 0.0im
    chain1.ϵ  != chain2.ϵ  && return 0.0im
    # 
    @unpack j0 = chain1
    # 
    return 1e6*sum(
        Osym(chain1,λ,τ,σs)'*
            Osym(chain2,λ,τ,σs)
                for λ in -j0:j0, τ ∈ -jω:jω)
end




JPs = [jp"0+"]#,jp"1+",jp"1-",jp"2+",jp"2-",jp"3+",jp"4+"]

# JP in JPs,
waveset = selectwaveset(jp"0+")
# waveset = selectwaveset(jp"1-")

function enhance(m::Matrix)
    d = Diagonal(1 ./ sqrt.(diag(m)))
    return d * m * d
end


# all
Nwaves = size(waveset,1)
Φij = zeros(Complex{Float64}, Nwaves,Nwaves)
for i in 1:Nwaves, j in 1:Nwaves
    ci = Chain(; waveset[i,:]...)
    cj = Chain(; waveset[j,:]...)
    Φij[i,j] = three_body_phase_space_integral(σs->crossedterm(σs, ci, cj), ms)
end

wavelabels = latexstring.(ξL.([Chain(; waveset[i,:]...) for i in 1:Nwaves]))
heatmap(1:Nwaves, 1:Nwaves, ifsmallgiveanan.(real.(enhance(Φij))),
    yticks=(1:Nwaves, wavelabels),
    xticks=(1:Nwaves, wavelabels), xrotation = 90,
    c=:diverging_bwg_20_95_c41_n256, clim=(-1,1),
    size=(Nwaves*40+60, Nwaves*40), bottom_margin=10mm,
    title=latexstring("\\mathrm{Sectors}\\,\\,J^{PC}:"*prod(" $(JP.j)^{$(JP.p)+}" for JP in JPs)))
#
