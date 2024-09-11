using LinearAlgebra

function newton_raphson() end
function solve() end

"""
    diagonality(A)

Determines whether matrix, `A` is strictly, diagonally dominant.
"""
function diagonality(A::Matrix)
    n, m    = size(A)
    diags   = zeros(eltype(A), n, m)
    for i ∈ 1:1:n, j ∈ 1:1:m
        if i != j
            diags[i, j] = A[i, j]
        end
    end
    return tr(A) >= sum(diags)
end

"""
    spectral_radius(A)

Finds the spectral radius of matrix, `A`.

# Notes
``ρ(\\mathbf{A}) = \\max|λ|``, where λ is the set of eigvalsvalues for `A` [burdenNumericalAnalysis2016]_.
"""
spectral_radius(A::Matrix) = maximum(abs.(eigvals(A)))

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
    At, is_symmetric = transpose(A), false
    i = 1; for ai ∈ A
        j = 1; for aj ∈ ai
            if aj == At[i, j]
                is_symmetric = true
            else
                is_symmetric = false
                return is_symmetric
            end
            j += 1
        end
        i += 1
    end
    return is_symmetric
end

"""
    positive_definite(A)

Determines whether matrix, `A` is positive definite.
"""
positive_definite(A::Matrix) = issymmetric(A) && all(eigvals(A) .>= 0.)

"""
    tridiagonality(A)

Determine whether matrix, `A` is tridiagonal.
"""
tridiagonality(A::Matrix) = sum(A - Tridiagonal(diag(A, -1), diag(A), diag(A, 1))) == 0.