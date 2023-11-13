using ThreeBodyDecay

const mηp = 0.958; const Γηp = 0.196e-3;
const mK0 = 1.4256; const ΓK0 = 0.270;
const mK = 0.493;


sqrt(0.3)+mK
ρηK(s) = quadgk(σ->sqrt(λ(s,σ,mK^2))*abs2(amp(σ,BreitWigner(mηp,Γηp))),0.68,(√s-mK)^2)[1]/s
ρηK_stable(s) = sqrt(λ(s,mηp^2,mK^2))/s

let
    n = ρηK(2.0)
    plot(e->ρηK(e^2),1.32,2)
    n = ρηK(2.0)
    plot!(e->sqrt(λ(e^2,mηp^2,mK^2))/e^2,mηp+mK,2)

IK0(s) = abs2(amp(s,BreitWigner(mK0,ΓK0))) * ρηK(s)
IK0_stable(s) = abs2(amp(s,BreitWigner(mK0,ΓK0)))*ρηK_stable(s)

let
    n = IK0(2.0^2)
    plot(e->IK0(e^2)/n,1.32,3, lab=L"\Gamma_{\eta} = 0.2\,\mathrm{MeV}")
    n = IK0_stable(2.0^2)
    plot!(e->IK0_stable(e^2)/n,mηp+mK,3, lab=L"\Gamma_\eta = 0")
    plot!(xlab=L"\sqrt{s}\equiv m_{\eta K}\,\mathrm{(GeV)}",
        ylab=L"|\mathrm{BW}_{K_{0}}|^2\Phi_{K\eta}(s)",
        title=L"K_{0}(1430)\to K\,\eta(958)\,\,S\!-\!\mathrm{wave}")
end
savefig("K02Ketap_v2.pdf")
