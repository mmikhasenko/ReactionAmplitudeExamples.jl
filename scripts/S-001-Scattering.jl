const mΣc = 2.4529
const mDbar = 1.86482

const mJψ = 3.092
const mp = 0.938

#
mJψ+mp
mΣc+mDbar

#
ρ1(s) = s < (mΣc+mDbar)^2 ? 0.0 : sqrt(λ(s,mΣc^2,mDbar^2))/s;
ρ2(s) = sqrt(λ(s,mJψ^2,mp^2))/s
#
M(s,γ) = 1/(-γ-1im*sqrt(s-(mΣc+mDbar)^2)/2)
pp(γ) = (mΣc+mDbar)^2-(2γ)^2

function plot_complex_plane_and_spectra(γ)
    plot(layout=grid(3,1,heights=(0.5,0.25,0.25)),size=(700,800), link=:x)
    plot!([4.25^2, 4.4^2], [0.0, 0.0], lab="", l=(2,:red), α=0.1, sp=1)
    plot!([(mDbar+mΣc)^2, 4.4^2], [0.0, 0.0], lab="", l=(2,:red), sp=1)
    plot!([(mDbar+mΣc)^2], [0.0], sp=1, m=(3,:red), lab="",
        ann=(4.33^2, 0.6, "[Dbar Sigmac]"))
    plot!([real(pp(γ))], [1000*imag(sqrt(pp(γ)))], sp=1,
        m=(7,(imag(pp(γ)) > 0 ? :orange : :purple)), ylim=(-5,5),
        xlab="mass square (GeV)", ylab="- half width (MeV)", lab="")
    #
    plot!(s->abs2(M(s+1e-6im,γ))*ρ1(s), 4.25^2, 4.4^2, yaxis=false, lab="",
        xlab="m(Dbar Sigma) (GeV)", l=(2,:black), sp=2)
    plot!(s->abs2(M(s+1e-6im,γ))*ρ2(s), 4.25^2, 4.4^2, yaxis=false, lab="",
        xlab="m(Jpsi p) (GeV)", l=(2,:black), sp=3)
end

# plot_complex_plane_and_spectra((200+25im)*1e-3)

@manipulate for γRe = LinRange(-150,350,100), γIm = LinRange(0,20,100)
    γ = (γRe + 1im*γIm)*1e-3
    plot_complex_plane_and_spectra(γ)
end

# @manipulate for Esq = LinRange(4.25^2,4.4^2,50), HG = LinRange(-5e-3,4e-3,50)
#     γ = sqrt((mΣc+mDbar)^2-(Esq+2im*HG*sqrt(Esq)))/2
#     plot_complex_plane_and_spectra(γ)
# end
