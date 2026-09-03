using Revise
using Pkg

includet("Anderson.jl")
using .Anderson

Pkg.add("Plots")
using Plots

using LinearAlgebra

Pkg.add("Krylov")
using Krylov

Pkg.add("MatrixDepot")
using MatrixDepot

function affine()
    n = 500
    md = mdopen("gravity", n, false)
    A = md.A
    b = md.b

    x0 = zeros(n)
    g = x -> x - (A * x - b)


    x_cg, stats_cg = gmres(A, b, x0, atol=1e-6)
    
    @show stats_cg.niter, size(stats_cg.residuals)
    
    x_a, stats_a  = anderson(x0, g, typemax(Int))

    @show stats_a.iterations, stats_a.residual_history[end]
end

affine()