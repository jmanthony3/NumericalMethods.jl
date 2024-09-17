module SystemOfEquations

export SystemOfEquation
export conjugate_gradient
export gaussian_elimination
export steepest_descent
# export solve

import ..NumericalMethods: solve, positive_definite, PositiveDefiniteError

using LinearAlgebra

"""
    SystemOfEquation(A, b, N, tol)

Solve a linear system of equations (SOE): ``\\mathbf{A}\\vec{x} = \\vec{b}``.
"""
struct SystemOfEquation{T<:Real}
    A   ::Matrix{T}
    b   ::Vector{T}
    N   ::Int64
    tol ::T
end

## Ch. 6.1 (p. 362)
"""
    gaussian_elimination(SOE)

Directly find the solution to linear `SOE` by *Gaussian Elimination with Back Substitution*.
"""
function gaussian_elimination(SOE::SystemOfEquation)::Vector{Float64}
    n = size(SOE.A)[1]
    Aug = zeros((n, n + 1))
    Aug[1:n,1:n]   .= SOE.A
    Aug[:,n+1]     .= SOE.b
    for i in 1:1:n
        E = Aug
        p = findfirst(Aug[i:n,i] .!= 0) + (i - 1)
        if i != p
            ei, ep = E[i, :], E[p, :]
            E[i, :] .= ep
            E[p, :] .= ei
        end
        for j in i+1:1:n
            mji = Aug[j,i] / Aug[i,i]
            E[j,:] .= E[j,:] - mji .* E[i,:]
        end
        Aug = E
    end
    if Aug[n,n] == 0
        error("No unique solution could be found.")
    end
    x = zeros(n)
    x[n] = Aug[n,n+1] / Aug[n,n]
    for i in n-1:-1:1
        aijxj = 0.
        for j in i+1:1:n
            aijxj += Aug[i,j] * x[j]
        end
        x[i] = (Aug[i,n+1] - aijxj) / Aug[i,i]
    end
    return x
    # return (k <= N && isassigned(approximations, k)) ? approximations[k] : NaN
end

## Ch. 7.6 (p. 487)
"""
    conjugate_gradient(SOE, x[, C])

Use initial guess vector, `x` and (if desired) pre-conditioning matrix, `C` to solve `SOE`.

*Is best suited for large, sparse matrices*.
If pre-conditioned, can solve in ``\\sqrt{n}`` iterations.
More computationally expensive than `gaussian_elimination` for smaller systems.
"""
function conjugate_gradient(SOE::SystemOfEquation,
        x   ::Vector{Float64},
        C   ::Union{Nothing, Bool, Matrix{Float64}}=nothing)::Vector{Float64}
    if positive_definite(SOE.A)
        r0              = SOE.b - (SOE.A * x)
        do_precondition = true
        v0              = r0
        Minv            = SOE.A
        C = if isnothing(C)
            # do_precondition = true
            # v0 = r0
            nothing
        elseif isa(C, Bool)
            if C
                do_precondition = true
                Minv = inv(diag(SOE.A) * transpose(diag(SOE.A)))
                v0 = Minv * r0
                diag(SOE.A)
            else
                do_precondition = false
                v0 = r0
                false
            end
        else
            do_precondition = false
            Minv = inv(C * transpose(C))
            v0 = (Minv * r0)
            C
        end
        k, approximations, errors = 1, [x], [SOE.tol*10]
        while last(errors) > SOE.tol && k <= SOE.N
            α = float(if do_precondition
                (transpose(r0) * r0)
            else
                ((transpose(r0) * Minv) * r0)
            end / ((transpose(v0) * SOE.A) * v0))
            x1 = x + α .* v0
            # push!(approximations, x1)
            push!(errors, norm(x1 - x))
            r1 = r0 - α .* (SOE.A * v0)
            s1 = float(if do_precondition
                (transpose(r1) * r1) / (transpose(r0) * r0)
            else
                ((transpose(r1) * Minv) * r1) / ((transpose(r0) * Minv) * r0)
            end)
            x, r0 = x1, r1
            v0 = (do_precondition ? r1 : (Minv * r1)) + (s1 .* v0)
            k += 1
        end
        return x
        # return (k <= N && isassigned(approximations, k)) ? approximations[k] : NaN
    else
        throw(PositiveDefiniteError)
    end
end

## Ch. 10.4 (p. 666)
"""
    steepest_descent(SOE, x)

Approximate solution to `SOE` with initial guess vector, `x`.
"""
function steepest_descent(SOE::SystemOfEquation, x::Vector{Float64})::Vector{Float64}
    k, approximations, errors = 1, [x], [SOE.tol*10]
    while last(errors) > SOE.tol && k <= SOE.N
        r = SOE.b - (SOE.A * x)
        α = float((transpose(r) * r) / ((transpose(r) * SOE.A) * r))
        x1 = x + α .* r
        # push!(approximations, x1)
        push!(errors, norm(x1 - x))
        x = x1; k += 1
    end
    return x
    # return (k <= N && isassigned(approximations, k)) ? approximations[k] : NaN
end

"""
    solve(SOE::SystemOfEquation[; method=:gaussian_elimination, x=nothing, C=nothing])

Solve ``\\mathbf{A}\\vec{x} = \\vec{b}`` according to `method` ∈ {`:gaussian_elimination` (default), `:steepest_descent`, `:conjugate_gradient`}.

Each `method` has an equivalent convenience function.
E.g. `solve(SOE; method=:gaussian_elimination)` ≡ `gaussian_elimination(SOE)`.
"""
function solve(SOE::SystemOfEquation;
    method  ::Symbol                                = :gaussian_elimination,
    x       ::Union{Nothing, Vector{Float64}}       = nothing,
    C       ::Union{Nothing, Bool, Matrix{Float64}} = nothing)::Vector{Float64}
    if size(SOE.A)[1] == size(vec(SOE.b))[1]
        return if method == :gaussian_elimination
            gaussian_elimination(SOE)
        elseif method == :steepest_descent
            steepest_descent(SOE, x)
        elseif method == :conjugate_gradient
            conjugate_gradient(SOE, x, C)
        end
    else
        throw(SystemOfEquationError(SOE.A, SOE.b))
    end
end

end
