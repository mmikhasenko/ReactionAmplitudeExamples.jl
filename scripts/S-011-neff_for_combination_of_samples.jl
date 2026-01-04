using Statistics
using Plots


N1 = 100
w1 = rand(N1)

N2 = 100
w2 = abs2.(randn(N2))

histogram(w1)
histogram(w2)

Neff(w) = sum(w)^2 / sum(w.^2)

Neff(w1)
Neff(w2)

Neff_scale(a) = Neff([w1..., (a .* w2)...])

Neff_scale(1.0)
Neff_scale(2.0)

plot(a->Neff_scale(a), 0.1, 1.5)
