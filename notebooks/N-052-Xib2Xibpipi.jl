### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 9929c400-4607-11ed-0cc7-db923ef1334a
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec("PyCall"),
        Pkg.PackageSpec("SymPy"),
        Pkg.PackageSpec("Parameters")
    ])
    # 
    using Parameters
    using SymPy
    # 
    import PyCall
    PyCall.pyimport_conda("sympy.physics.wigner", "sympy")
    PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
    # 
    import_from(sympy.physics.wigner)
    # import_from(sympy.physics.quantum.cg, )
    import_from(sympy.physics.quantum.spin, (:WignerD, :CG), typ=:Any)
end

# ╔═╡ 072ab567-6b5b-4e90-9127-1debf2a9cffa
md"""
## Transition $\text{baryon} \to \text{baryon}\,0^-\,0^-$
via intermediate $1/2$ and $3/2$
"""

# ╔═╡ ee76c86b-0d6c-423d-9a5c-870821bb68c8
# 2022-10-07 Misha Mikhaseko

# ╔═╡ 0c1ce318-7106-41fc-94a1-af5cd0e1e738
begin
    @with_kw struct Chain{T}
        j::T
        l::T
    end
    # 
    function Chain(jp::String)
        two_j = jp[1] - '0'
        iseven(two_j) && error("two_j is not odd! Expect baryon")
        j = two_j / Sym(2)
        filter_l = jp[end] == '-' ? iseven : isodd
        #
        possible_l = div.((two_j .+ [-1, 1]), 2)
        selected_l = filter(filter_l, possible_l)
        l = selected_l[1] |> Sym
        Chain(; j, l)
    end
end

# ╔═╡ bfdb72ec-c298-4584-a9fb-abbabe949c3b
@syms θ::positive => "\\theta"

# ╔═╡ 8a4791c4-56a1-4404-b9aa-5dd67151deeb
dl(j0, l, λ, ν; ξ::Chain) =
    sqrt((2l + Sym(1)) * (2ξ.l + Sym(1)) / (2j0 + Sym(1))) *
    CG(l, 0, ξ.j, ν, j0, ν) * WignerD(ξ.j, ν, λ, 0, θ, 0) * CG(ξ.l, 0, 1 / Sym(2), λ, ξ.j, λ)

# ╔═╡ 90156e3d-4c68-4069-99a5-941a23682342
let # call the function for symbols
    @syms λ ν j l L j0 => "j_0"
    I = dl(j0, L, λ, ν; ξ=Chain(; j, l))
    # 
    Markdown.parse(
        "```math\n\\begin{align}\n \\mathcal{A}_{" *
        sympy.latex(ν) * "," * sympy.latex(λ) *
        "} &=" *
        sympy.latex(I) *
        "\\end{align}\n```")
end

# ╔═╡ 550a063a-34e8-4ab8-9659-9ac126a08b5d
function I(j0, l; ξ::Chain)
    sum(abs2(dl(j0, l, two_ν / Sym(2), two_λ / Sym(2); ξ).doit())
        for two_ν in -2j0:2:2j0,
        two_λ in -2ξ.j:2:2ξ.j)
end

# ╔═╡ 9120bb6a-d687-431a-aa98-5d3b60fefda7
md"""
## Validation with $\Omega_b^- \to \Xi_c^+ K^- \pi^-$
"""

# ╔═╡ 7ecd9d76-583c-4d08-9332-3a80307227cc
md"""
Cross check is done with the $\Omega_b^- \to \Xi_c^+ K^- \pi^-$,
where we know (see [LHCb-PAPER-2021-012] (https://lhcbproject.web.cern.ch/Publications/LHCbProjectPublic/LHCb-PAPER-2021-012.html)) that the angular discributions for the decays of $\Omega_c^{**0}$ decaying to $\Xi_c^+ K^-$ are non-trivial.

The decays are:
 - b-decay $\Omega_b^{-} 1/2^\pm (\text{weak})$
 - two-body $\Xi_c^{+} (1/2^+)$ and $K^- (0^-)$ $\Rightarrow 1/2^{-}$ in S-wave
 - The $\Omega_c^{**0}$ states in P-multiplet have quantum numbers $1/2^-$, $3/2^-$, and $5/2^-$

| $\Omega_c^{**0}$ $J^P$ |  $\Xi_c^+ K^-$ $L$-wave | $I(\cos\theta)$|
|---|---|---|
| $1/2^-$ | $S$ | $c_0$ |
| $3/2^-$ | $D$ | $c_1\cos(\theta)^2 + c_0$ |
| $5/2^-$ | $D$ | $c_2\cos(\theta)^4 + c_1\cos(\theta)^2  + c_0$ |

"""

# ╔═╡ b76faacb-5756-47fe-be4f-7a6c76c488b3
md"""
The result of calculations in this notebook follows:
"""

# ╔═╡ 87c3caed-8085-4de7-a0a1-c65db729229d
Ob2XcKpi_intensity = [
    I(1 / Sym(2), 0; ξ=Chain("1/2-")),
    I(1 / Sym(2), 1; ξ=Chain("3/2-")),
    I(1 / Sym(2), 2; ξ=Chain("5/2-"))
] .|> simplify

# ╔═╡ c48dfb50-b64b-4ef3-9545-d1f301f25ab6
Markdown.parse("""
```math
\\begin{align}
""" *
               "&1/2^-: &&" * sympy.latex(Ob2XcKpi_intensity[1]) * "\\\\" *
               "&3/2^-: &&" * sympy.latex(Ob2XcKpi_intensity[2]) * "\\\\" *
               "&5/2^-: &&" * sympy.latex(Ob2XcKpi_intensity[3]) * "\\\\" *
               """
               \\end{align}
               ```
               """)

# ╔═╡ ecb054fb-438c-456e-af6c-9f65ac3f4616
md"""
## The study case: $\Xi_b^{**0,-} \to \Xi_b^{0,-} \pi^+ \pi^-$:
"""

# ╔═╡ 8842f56d-5d3a-4aa0-890f-63657e1cf07e
begin
    Ξc′ = Chain("1/2+")
    Ξcˣ = Chain("3/2+")
end;

# ╔═╡ ebc8e328-1ce5-4a56-81fc-29f37ca17d28
via_Ξc′ = [
    I(1 / Sym(2), 0; ξ=Ξc′),
    I(3 / Sym(2), 2; ξ=Ξc′)
] .|> simplify

# ╔═╡ 89de2412-ee68-44ce-b209-85677d96e806
via_Ξcˣ = [
    I(3 / Sym(2), 0; ξ=Ξcˣ),
    I(1 / Sym(2), 2; ξ=Ξcˣ)
] .|> simplify

# ╔═╡ Cell order:
# ╟─072ab567-6b5b-4e90-9127-1debf2a9cffa
# ╠═ee76c86b-0d6c-423d-9a5c-870821bb68c8
# ╠═9929c400-4607-11ed-0cc7-db923ef1334a
# ╠═0c1ce318-7106-41fc-94a1-af5cd0e1e738
# ╠═bfdb72ec-c298-4584-a9fb-abbabe949c3b
# ╠═8a4791c4-56a1-4404-b9aa-5dd67151deeb
# ╟─90156e3d-4c68-4069-99a5-941a23682342
# ╠═550a063a-34e8-4ab8-9659-9ac126a08b5d
# ╟─9120bb6a-d687-431a-aa98-5d3b60fefda7
# ╟─7ecd9d76-583c-4d08-9332-3a80307227cc
# ╟─b76faacb-5756-47fe-be4f-7a6c76c488b3
# ╠═87c3caed-8085-4de7-a0a1-c65db729229d
# ╟─c48dfb50-b64b-4ef3-9545-d1f301f25ab6
# ╟─ecb054fb-438c-456e-af6c-9f65ac3f4616
# ╠═8842f56d-5d3a-4aa0-890f-63657e1cf07e
# ╠═ebc8e328-1ce5-4a56-81fc-29f37ca17d28
# ╠═89de2412-ee68-44ce-b209-85677d96e806
