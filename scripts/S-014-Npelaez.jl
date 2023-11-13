module GKPY

# using Plots
using Parameters
using QuadGK

# masses from the paper
const mπ = 0.13957;
const mρ = 0.7736;
const Γρ = 0.146;
const mK = 0.496;

export mπ

# fit parameters (constrained)
const Ppars = (B0 = 1.043, B1 = 0.19, λ1 = 1.39, λ2 = -1.70, ϵ1 = 0.00, ϵ2 = 0.07, e0 = 1.05)
const Gpars = (B0 = 1.09e5, B1 = 1.41e5, λ = 0.051e5)

export δ1, δ3
export cotδ1, cotδ3

# maps
k(s) = sqrt(s/4-mπ^2)
w(s; s0 = 1.45^2) = (sqrt(s) - sqrt(s0-s)) / (sqrt(s) + sqrt(s0-s))

# phase shifts
# P-wave
function cotδ1_less_1050(s; pars)
    s < 4mπ^2 && return 0.0
    @unpack e0, B0, B1 = pars
    return sqrt(s)/(2*k(s)^3)*(mρ^2-s)*(2mπ^3/(mρ^2*sqrt(s)) + B0 + B1*w(s; s0 = e0^2))
end
function δ1_more_1050(s; pars)
    s < 4mπ^2 && return 0.0
    @unpack  λ1, λ2, ϵ1, ϵ2 = pars
    λ0 = acot(cotδ1_less_1050(4mK^2; pars=pars))
    return λ0 + λ1*(sqrt(s)/(2mK)-1) + λ2*(sqrt(s)/(2mK)-1)^2
end
function _δ1(s; pars = Ppars)
    @unpack e0 = pars
    if s < e0^2 
        v = acot(cotδ1_less_1050(s; pars=pars))
        return (v<0) ? v+π : v
    end
    return s < 1.4^2 ? δ1_more_1050(s; pars=pars) : δ1_more_1050(1.4^2; pars=pars)
end

function δ1(s; pars = Ppars)
    v = _δ1(s; pars = Ppars)
    s<0.8^2 ? v : (v>0 ? v : v+π)
end


cotδ1(s; pars = Ppars) = cot(δ1(s; pars=pars))

end


using Main.GKPY
using Plots
using QuadGK


theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right)


plot(e->δ1(e^2), 2mπ, 3.4)


t(s) = sin(δ1(s))*exp(1im*δ1(s))
let
    plot()
    plot!(e->real(t(e^2)), 2mπ, 1.4)
    plot!(e->imag(t(e^2)), 2mπ, 1.4)
end

Omnes(s) = exp(s*quadgk(s′->δ1(s′)/s′/(s′-s-1e-6im), 4mπ^2, Inf)[1]/π)

let
    plot()
    plot!(e->δ1(e^2), 2mπ, 3.4)
    plot!(e->imag(log(Omnes(e^2))), 2mπ, 3.4)
end
using LaTeXStrings
let
    plot(ylab=L"t(e^2) / \Omega(e^2)", xlab=L"e\,\,[\mathrm{GeV}]", leg=:left)
    plot!(e->real(t(e^2)/Omnes(e^2)), 2mπ, 1.4, lab="Re")
    plot!(e->imag(t(e^2)/Omnes(e^2)), 2mπ, 1.4, lab="Im")
end