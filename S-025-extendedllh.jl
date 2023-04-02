using Symbolics
using LinearAlgebra

N = 10
@variables xs[1:10]
xs = collect(xs)

function extendednegativeloglikelihood(xs, ϕs_all,M)
    μ = xs' * M * xs
    logsum = sum(ϕs_all) do ϕs
        log(abs2(xs' * ϕs))
    end
    return -logsum + μ
end

ellh = extendednegativeloglikelihood(xs, [rand(10) for i in 1:3], rand(10,10))

ellh

grad = Symbolics.gradient(ellh, xs)
hes = Symbolics.jacobian(grad, xs)
hes_sparse = Symbolics.sparsejacobian(grad, xs)
hes_sparsity = Symbolics.jacobian_sparsity(grad, xs)
