using Revise
using Pkg

includet("Anderson.jl")
using .Anderson

Pkg.add("Plots")
using Plots

function bratu(u; λ = 1.0)
    n = length(u); h = 1 / (n + 1)
    v = similar(u)
    for i in 1:n
        ul = i == 1 ? 0.0 : u[i-1]
        ur = i == n ? 0.0 : u[i+1]
        v[i] = 0.5 * (ul + ur + h^2 * λ * exp(u[i]))
    end
    v
end

g = x -> bratu(x; λ = 1.0)
x0 = zeros(500)

x, stats = anderson(x0, g, 60, max_iter=5000)

@show stats.iterations, stats.residual_history[end]
