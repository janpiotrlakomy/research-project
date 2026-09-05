module Anderson

export anderson, scalar_newton


import LinearAlgebra: norm

function anderson(x0::AbstractVector{T}, g, m::Integer; β=1.0::AbstractFloat, tol=1e-6::AbstractFloat, max_iter=100::Integer) where T <: AbstractFloat
    f = x -> g(x) - x

    x_prev = copy(x0)
    f_prev = f(x_prev)

    x_curr = x_prev + β * f_prev
    f_curr = f(x_curr)

    X = Matrix{T}(undef, length(x0), 0)
    F = Matrix{T}(undef, length(x0), 0)

    history = [copy(x0), copy(x_curr)]
    residual_history = [norm(f(x0)), norm(f_curr)]

    for j in 1:max_iter
        if j > m 
            X = X[:, 2:end]
            F = F[:, 2:end]
        end

        dx = x_curr - x_prev
        df = f_curr - f_prev

        X = hcat(X, dx)
        F = hcat(F, df)

        # solve f_j = F_j * γ for γ
        γ = F \ f_curr



        x_prev = x_curr
        f_prev = f_curr

        x_curr = (x_prev - X * γ) + β * (f_prev - F * γ)
        f_curr = f(x_curr)

        push!(residual_history, norm(f_curr))
        push!(history, x_curr)
        if (norm(f_curr) < tol)
            return (x = x_curr,stats=( residual = norm(f_curr), iterations = j,  history = history, residual_history = residual_history))
        end
    end

    return (x = x_curr,stats=( residual = norm(f_curr), iterations = max_iter, history = history, residual_history = residual_history))
    

end # function anderson


function scalar_newton(x0::AbstractFloat, g, dg, tol=1e-6::AbstractFloat, max_iter=100::Integer)
    f = x -> g(x) - x
    df = x -> dg(x) - 1.0

    history = [x0]
    residual_history = [norm(f(x0))]
    x = x0
    for i in 1:max_iter
        if norm(f(x)) < tol
            return (x = x, stats=(residual = norm(f(x)), iterations = i, history = history, residual_history = residual_history))
        end

        x = x - f(x) / df(x)
        push!(history, x)
        push!(residual_history, norm(f(x)))
    end
    return (x = x, stats=(residual = norm(f(x)), iterations = max_iter, history = history, residual_history = residual_history))
end # function scalar_newton

end # module Anderson 