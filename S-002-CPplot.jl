const mΣc = 2.4529
const mDbar = 1.86482

const mJψ = 3.092
const mp = 0.938

pyplot()
let
    plot(xlab=L"\mathrm{Re}\,s\,\,(\mathrm{GeV}^2)",
         ylab=L"\mathrm{Im}\, s\,\,(\mathrm{GeV}^2)",size=(500,350))
    plot!([(mJψ+mp)^2, 4.4^2], [0.0, 0.0], lab="", l=(2,:red), sp=1,
        ylim=(-1,1), xlim = ((mJψ+mp)^2-1, 4.4^2),
            ann=((mJψ+mp)^2, 0.1, L"[\,J/\psi\,p\,]"))
    plot!([(mJψ+mp)^2], [0.0], lab="", m=(8,:red), sp=1)
    #
    plot!([4.312^2], [-0.5], sp=1, m=(15, :purple), lab="",
        ann=(4.3^2, -0.5, text(L"P_c(4312)",:right)))
end
savefig("figs/Jpsi_pole.pdf")
