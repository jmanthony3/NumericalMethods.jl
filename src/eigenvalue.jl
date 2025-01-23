module Eigenvalues

export Eigenvalue
export power_method
export inverse_power_method
export qr_algorithm

import LUSE_ENGR701_704_NumericalMethods: symmetry, SymmetricError, tridiagonality

using LinearAlgebra

# Ch. 9 (p. 569)
"""
    Eigenvalue(A, N, tol)

Find the characteristic (eigen) values of matrix, `A`.
Typically thought as roots of polynomial from determinant.
"""
struct Eigenvalue{T<:Real}
    A   ::Matrix{T}
    N   ::Int64
    tol ::T
end

## 9.3 (p. 585)
"""
    power_method(EV, x0)

Approximate the dominant eigenvalue of matrix, `EV.A` given some non-zero, initial guess vector, `x0`.
"""
function power_method(EV::Eigenvalue,
        x0::Vector{Float64})::Float64
    μs  = [norm(x0, Inf)]
    x   = x0 ./ last(μs)
    k, eigenvectors, errors = 1, [x], [EV.tol*10]
    while last(errors) > EV.tol && k <= EV.N
        yp, y = 0., EV.A * x
        for yi in y
            if abs(yi) == norm(y, Inf)
                yp = yi
            end
        end
        eigenvector = y ./ yp
        push!(μs, yp)
        # push!(eigenvectors, eigenvector)
        push!(errors, norm(x - eigenvector, Inf))
        x = eigenvector; k += 1
    end
    return last(μs)
end

"""
    inverse_power_method(EV, x0, q)

Approximate eigenvalue closest to target, `q` of matrix, `EV.A` given some non-zero, initial guess vector, `x0`.
*Supposed to converge faster than `power_method` [burdenNumericalAnalysis2016]_.*
"""
function inverse_power_method(EV::Eigenvalue,
        x0::Vector{Float64}, q::Float64)::Float64
    A = inv(EV.A - q .* I(size(EV.A)[1]))
    x = x0
    μs = [1/norm(x, Inf) + q]
    k, eigenvectors, errors = 1, [x], [EV.tol*10]
    while last(errors) > EV.tol && k <= EV.N
        yp, y = 0., A * x
        for yi in y
            if abs(yi) == norm(y, Inf)
                yp = float(yi)
            end
        end
        eigenvector = y ./ yp
        push!(μs, 1/yp + q)
        # push!(eigenvectors, eigenvector)
        push!(errors, norm(x - eigenvector, Inf))
        x = eigenvector; k += 1
    end
    return last(μs)
end

## 9.5 (p. 610)
"""
    qr_algorithm(EV)

Systematically find all the eigenvalues of matrix, `EV.A`.

**Matrix must be symmetric and tridiagonal!**

This method is preferred over `power_method` and `inverse_power_method` by keeping round-off error to a minimum [burdenNumericalAnalysis2016]_.
Refer to this [example](https://www.youtube.com/watch?v=FAnNBw7d0vg) for an explanation and demonstration.
"""
function qr_algorithm(EV::Eigenvalue)::Vector{Float64}
    if issymmetric(EV.A) && tridiagonality(EV.A)
        A = EV.A
        m, n = size(A)
        k, eigenvectors, errors = 1, [diag(A)], [EV.tol*10]
        while last(errors) > EV.tol && k <= EV.N
            Q = zeros((m, n))
            R = zeros((m, n))
            QI = []
            for j in 1:1:n
                ai = zeros(m)
                for i in 1:1:m
                    ai[i] = A[i,j]
                end
                ai_perp = zeros(m)
                for i in 1:1:j-1
                    R[i,j] = dot(ai, QI[i])
                    ai_perp .= ai_perp + R[i,j] .* QI[i]
                end
                ai = ai - ai_perp
                R[j,j] = √(sum(ai .^ 2.))
                qi = ai ./ R[j,j]
                push!(QI, qi)
                for (i, q) in enumerate(qi)
                    Q[i,j] = q
                end
            end
            A = R * Q
            push!(eigenvectors, diag(A))
            push!(errors, 2\sum([norm(diag(A, -1)), norm(diag(A, 1))]))
            k += 1
        end
        return last(eigenvectors)
    else
        throw(SymmetricError)
    end
end

end
