using Revise
using Pkg

includet("Anderson.jl")
using .Anderson

Pkg.add("Plots")
using Plots

function anderson_vs_newton(x0::AbstractFloat, g, g_str::AbstractString, dg, m=10::Integer, β = 1.0::AbstractFloat, tol=1e-6::AbstractFloat, max_iter=100::Integer)
    res = scalar_newton(x0, g, dg)
    res_anderson = anderson([x0], x -> [g(x[1])], 10) 
    anderson_ys = only.(res_anderson.residual_history)

    newton_ys = only.(res.residual_history)

    plot(1:length(anderson_ys),
        anderson_ys .+ ones(length(anderson_ys)) .* 1e-10,
        yscale=:log10,
        lc=:red, 
        lw=1, 
        alpha=0.6, 
        label = "metoda Andersona"
    )

    plot!(1:length(newton_ys), 
        newton_ys.+ ones(length(newton_ys)) .* 1e-10,yscale=:log10,
        lc=:green, 
        lw=1, 
        alpha=0.6, 
        label = "metoda Newtona"
    )

    title!("Porównanie zbieżności metod Andersona i Newtona dla funkcji g(x) = $g_str", titlefontsize=8)

end

 f = x -> x^3 - 4*x^2 + 3*x + 2
    g = x -> x - f(x)
    dg = x -> 1 - 3*x^2 + 8*x - 3

    x0 = 0.0


anderson_vs_newton(x0, g, "-x^3 + 4*x^2 - 2*x - 2", dg)

g = x -> cos(x)
dg = x -> -sin(x)
x0 = 0.0

anderson_vs_newton(x0, g, "cos(x)", dg)