using LinearAlgebra

function newton_raphson() end
function solve() end

function issquare(A::Matrix)
    m, n = size(A)
    if m != n
        throw(SquareMatrixError((m, n)))
    end
    return true
end

"""
    diagonality(A)

Determines whether matrix, `A` is strictly, diagonally dominant.
"""
function diagonality(A::Matrix)
    if issquare(A)
        m, n    = size(A)
        diags   = zeros(eltype(A), m, n)
        for i ∈ 1:1:m, j ∈ 1:1:n
            if i != j
                diags[i, j] = A[i, j]
            end
        end
        return tr(A) >= sum(diags)
    end
end

"""
    spectral_radius(A)

Finds the spectral radius of matrix, `A`.

# Notes
``ρ(\\mathbf{A}) = \\max|λ|``, where λ is the set of eigvalsvalues for `A` [burdenNumericalAnalysis2016]_.
"""
spectral_radius(A::Matrix, N::Int64=100, tol::Real=10^-6)   = maximum(abs.(if issquare(A) && issymmetric(A) && tridiagonality(A)
    qr_algorithm(Eigenvalue(A, N, tol))
else
    eigvals(A)
end))
# spectral_radius(A::Matrix)                                  = maximum(abs.(eigvals(A)))

"""
    condition_number(A)

Find the condition number of matrix, `A`.

# Notes
Definition [burdenNumericalAnalysis2016]_:
The condition number of the non-singular matrix, ``\\mathbf{A}`` relative to a norm, ||⋅|| is

```math
    K(\\mathbf{A}) = ||\\mathbf{A}|| ⋅ ||\\mathbf{A}^{-1}||
```

A matrix is well-conditioned if ``K(\\mathbf{A})`` is close to 1 and is ill-conditioned if significantly greater than 1.
"""
condition_number(A::Matrix) = norm(A) * norm(inv(A))

"""
    symmetry(A)

Determines whether matrix, `A` is symmetric.
"""
function symmetry(A::Matrix)
    if issquare(A)
        At, is_symmetric = transpose(A), false
        i = 1; for ai ∈ A
            j = 1; for aj ∈ ai
                if aj == At[i, j]
                    is_symmetric = true
                else
                    return false
                end
                j += 1
            end
            i += 1
        end
        return is_symmetric
    end
end

"""
    positive_definite(A)

Determines whether matrix, `A` is positive definite.
"""
positive_definite(A::Matrix, N::Int64=100, tol::Real=10^-6) = issymmetric(A) && all(if tridiagonality(A)
    qr_algorithm(Eigenvalue(A, N, tol))
else
    eigvals(A)
end .>= 0.)
# positive_definite(A::Matrix)                                = issymmetric(A) && all(eigvals(A) .>= 0.)

"""
    tridiagonality(A)

Determine whether matrix, `A` is tridiagonal.
"""
tridiagonality(A::Matrix) = sum(A - Tridiagonal(diag(A, -1), diag(A), diag(A, 1))) == 0.