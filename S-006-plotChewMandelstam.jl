using Plots
using LaTeXStrings
pyplot()

using ThreeBodyDecay

font = Plots.font("sans-serif", 15)
fontt = Plots.font("sans-serif", 18)
pyplot(guidefont=font, tickfont=font, legendfont=font, titlefont = fontt)

import PyCall
import PyPlot: matplotlib
# font
PyCall.PyDict(matplotlib["rcParams"])["font.family"] = "serif"
PyCall.PyDict(matplotlib["rcParams"])["text.usetex"] = true
PyCall.PyDict(matplotlib["rcParams"])["text.latex.unicode"] = true

ChewMandestam(0.77^2+1e-5im, 0.15^2, 0.17^2)
iRho(0.77^2+1e-5im, 0.15^2, 0.17^2)

let sv = LinRange(-0.1, 1.4, 200), m1=0.14, m2=0.547
    cal1 = ChewMandestam.(sv .+ 1e-4im, m1^2, m2^2)
    cal2 = iRho.(sv .+ 1e-4im, m1^2, m2^2)
    plot(layout=grid(1,2), size=(900,350), ylim=(-1,1), grid=false, frame=:frame)
    plot!(sv,  [real.(cal2) imag.(cal2)], sp=1,
        ylab=L"16\pi\,\,i\rho(s)", xlab=L"s\,\,(\mathrm{GeV}^2)",
        lab=[L"\mathrm{Real\,\,part}" L"\mathrm{Imag\,\,part}"],
        ls=[:dash :solid], lc=[:black :red], lw=1.5)#, lc=[palette(:auto)[1] palette(:auto)[3]])
    plot!(sv, [real.(cal1) imag.(cal1)], sp=2,
        ylab=L"16\pi\,\,\Sigma(s)", xlab=L"s\,\,(\mathrm{GeV}^2)",
        lab=[L"\mathrm{Real\,\,part}" L"\mathrm{Imag\,\,part}"],
        ls=[:dash :solid], lc=[:black :red], leg=:bottomright, lw=1.5)#, lc=[palette(:auto)[1] palette(:auto)[3]])
    for sp=1:2
        vline!(sp=sp, [(m1+m2)^2], lab="", l=(1,0.5,:grey),
            ann=((m1+m2)^2, 0.5, text(L"s_{\mathrm{thr}}",:right,rotation=90)))
        hline!(sp=sp, [0], lab="", l=(1,0.5,:grey))
    end
    plot!()
end

# gr()
# let m1=0.14, m2=0.547
#     sx = LinRange(-0.1, 1.4, 200)
#     sy = LinRange(-1,1,100)
#     #
#     iRho2(s) = iRho(s, m1^2,m2^2)*s
#     cal1 = iRho2.(transpose(sx) .+ 1im .* sy)
#     heatmap(sx, sy, real.(cal1))
# end

savefig(joinpath("plot","ChewMandelstam.pdf"))
