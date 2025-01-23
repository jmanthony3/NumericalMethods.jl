module MultiVariableIterations

export MultiVariableIteration
export find_omega
export jacobian_form
# export solve
export gauss_seidel
export jacobi
# export newton_raphson
export successive_relaxation

import ..LUSE_ENGR701_704_NumericalMethods: issquare, SystemOfEquationError, solve, newton_raphson

using DataFrames # KEY: [JSSv107i04]
using LinearAlgebra
# KEY: [10.1145/3511528.3511535]
using Symbolics: Num, @variables, Differential, simplify, expand_derivatives, build_function, value

# Ch. 7 (p. 437)
"""
    MultiVariableIteration(A, x, b, N, tol)

Given the system of equations, ``\\mathbf{A}\\vec{x} = \\vec{b}``, attempt to solve ``\\vec{x} = \\mathbf{A}^{-1}\\vec{b}`` in so many iterations, `N` within tolerance, `tol`.
*Ideal for large, sparse systems.*
"""
struct MultiVariableIteration{T<:Real}
    A   ::Union{Matrix{T}, Vector{Function}}
    x   ::Vector{T}
    b   ::Vector{T}
    N   ::Int64
    tol ::T
end

__find_xk(T, x, c) = T * x + c

"""
    find_omega(MVI[, omega=0.])

Find the optimum relaxation parameter, `omega` for `MVI.A` and `MVI.x`, if possible.

# Notes
If 0 < `omega` < 2, then method will converge regardless of choice for `MVI.x`.
If matrix, `MVI.A` is tridiagonal, then the spectral radius of Gauss-Seidel's T-matrix, ``\\mathbf{T}\\_{g}`` is used to calculate ``ω := 2 / (1 + \\sqrt{(1 - `spectral_radius`(\\mathbf{T}\\_{g})})`` [burdenNumericalAnalysis2016]_.
"""
function find_omega(MVI::MultiVariableIteration, omega::Float64=0.)::Float64
    n, m = size(MVI.A)
    xn, xt = MVI.x, transpose(MVI.x)
    y = (xt * MVI.A) * xn
    theorem_6_22 = issymmetric(MVI.A) && y > 0.
    i, theorem_6_25 = 2, true; while i <= n && theorem_6_25 == true
        theorem_6_25 = det(MVI.A[begin:i, begin:i]) > 0; i += 1
    end
    return if theorem_6_22 || theorem_6_25
        if tridiagonality(MVI.A)
            D   = diagm(diag(MVI.A))
            L   = diagm(-1 => diag(MVI.A, -1))
            U   = diagm(1  => diag(MVI.A, 1))
            Tg  = inv(D - L) * U
            2. / (1. + √(1 - spectral_radius(Tg)))
        else
            omega
        end
    else
        omega
    end
end

function jacobian_form(g, variables::Tuple{Vararg{Num}})
    n = length(g)
    jacMatrix = Matrix{Function}(undef, n, n)
    for i ∈ 1:1:n, j ∈ 1:1:n
        t = variables[j]
        Dt = Differential(t)
        jacMatrix[i, j] = build_function(simplify(expand_derivatives(Dt(g[i](variables...))); expand=true), variables..., expression=Val{false})
    end
    return jacMatrix # * variables
end

function showjacobian(J::Matrix{Function}, x::Union{Tuple{Vararg{Num}}, Vector{Float64}})
    n, m = size(J)
    jacMatrix = Matrix{Num}(undef, n, m)
    for i ∈ 1:1:n, j ∈ 1:1:m
        jacMatrix[i, j] = J[i, j](x...)
    end
    return jacMatrix
end

"""
    solve(MVI::MultiVariableIteration[; method=:jacobi, omega=0., variables, jacobian])

Solve ``\\vec{x} = \\mathbf{A}^{-1}\\vec{b}`` according to `method` ∈ {`:jacobi` (default), `:gauss_seidel`, `:successive_relaxation`, `:newton_raphson`}.

Each `method` has an equivalent convenience function.
E.g. `solve(MVI; method=:jacobi)` ≡ `jacobi(MVI)`.
"""
function solve(MVI::MultiVariableIteration;
        method      ::Symbol                            = :jacobi,
        omega       ::Float64                           = 0.,
        variables   ::Union{Nothing, Tuple{Vararg{Num}}}= nothing,
        jacobian    ::Union{Nothing, Function}          = nothing)::Vector{Float64}
    if size(MVI.A)[1] == size(vec(MVI.x))[1] == size(vec(MVI.b))[1]
        x = MVI.x
        if method == :newton_raphson
            J(x) = if isnothing(jacobian)
                convert.(Float64, value.(showjacobian(jacobian_form(MVI.A, variables), x)))
            else
                jacobian(x)
            end
        elseif method == :successive_relaxation
            if omega == 0.
                ω = find_omega(MVI)
            # elseif isinstance(omega, (int, float)) && omega > 0.
            #     w = self.find_omega(omega=omega)
            #     logging.info(f"omega = {omega} given. Which is not optimum: {w}")
            #     w = omega
            else
                ω = omega
            end
        end
        k, n, approximations, errors = 1, length(x), [x], [MVI.tol*10]
        while last(errors) > MVI.tol && k <= MVI.N
            xk = zeros(eltype(MVI.x), n)
            for i ∈ 1:1:n
                if method == :newton_raphson
                    g = [f(x...) for f ∈ MVI.A]
                    xk = x + (J(x) \ -g)
                else
                    xk[i] = (if method ∈ (:gauss_seidel, :successive_relaxation)
                        y1 = sum([MVI.A[i, j] * xk[j] for j ∈ 1:1:i])
                        y2 = 0.; for j ∈ i+1:1:n
                            y2 += MVI.A[i, j] * x[j]
                        end
                        -y1 - y2
                    elseif method == :jacobi
                        -sum([j != i ? MVI.A[i, j] * x[j] : 0. for j ∈ 1:1:n])
                    end + MVI.b[i]) / MVI.A[i, i]
                end
            end
            if method == :successive_relaxation
                xk = [(1. - ω) * x[i] + ω * xk[i] for i ∈ 1:1:n]
            end
            # push!(approximations,   xk)
            push!(errors,           norm(xk - x))
            x = xk; k += 1 # iterate to k + 1
        end
        return x
        # return (k <= N && isassigned(approximations, k)) ? approximations[k] : NaN
    else
        throw(SystemOfEquationError(MVI.A, MVI.x))
    end
end

"""
    gauss_seidel(MVI)

Solve ``\\vec{x} = \\mathbf{A}^{-1}\\vec{b}`` via the **Gauss-Seidel Method**.

# Notes
This method improves on `jacobi()` by using the most recently calculated entries in the approximation vector, `x` at the end of each iteration.
The core algorithm by which method marches through iterations:
```math
    \\vec{x}^{(k)} = \\bigl( (\\mathbf{D} - \\mathbf{L})^{-1} * \\mathbf{U} \\bigr) ⋅
        \\vec{x}^{(k - 1)} + \\bigl( (\\mathbf{D} - \\mathbf{L})^{-1} \\bigr) ⋅ \\vec{b}
```
"""
gauss_seidel(MVI::MultiVariableIteration) = solve(MVI; method=:gauss_seidel)

"""
    jacobi(MVI)

Solve ``\\vec{x} = \\mathbf{A}^{-1}\\vec{b}`` via the **Jacobi Method** to find ``\\vec{x}``.

# Notes
The core algorithm by which method marches through iterations:
```math
    \\vec{x}^{(k)} = \\bigl( \\mathbf{D}^{-1} * (\\mathbf{L} + \\mathbf{U}) \\bigr) ⋅
        \\vec{x}^{(k - 1)} + ( \\mathbf{D}^{-1} ) ⋅ \\vec{b}
```
"""
jacobi(MVI::MultiVariableIteration) = solve(MVI; method=:jacobi)

"""
    newton_raphson(MVI::MultiVariableIteration, variables::Tuple{Vararg{Num}}[; jacobian=nothing])

Solve non-linear systems of equations, ``\\vec{x} = \\mathbf{A}^{-1}\\vec{b}`` via the **Newton-Raphson Method**.

Here, `MVI.A` should be a vector of functions wherein each variable is represented.
**Method will go faster if `jacobian` is pre-defined.**
_Otherwise, the Jacobian matrix of `MVI.A` will be internally constructed._

# Examples
```jldoctest; output=false
f1(x1, x2, x3)  = 3x1 - cos(x2*x3) - 0.5
f2(x1, x2, x3)  = x1^2 - 81(x2 + 0.1)^2 + sin(x3) + 1.06
f3(x1, x2, x3)  = exp(-x1*x2) + 20x3 + 3\\(10π - 3)
A               = [f1, f2, f3]
b               = zeros(length(A))
x0              = [0.1, 0.1, -0.1]
tol             = 1e-9
MVI             = MultiVariableIteration(A, x0, b, 5, tol)
using Symbolics
@variables x1, x2, x3
newton_raphson(MVI, (x1, x2, x3))

# output

3-element Vector{Float64}:
  0.5
 -1.176161909134556e-19
 -0.5235987755982989
```
"""
newton_raphson(MVI::MultiVariableIteration, variables::Tuple{Vararg{Num}};
    jacobian::Union{Nothing, Function}=nothing) = solve(MVI;
method=:newton_raphson, variables=variables, jacobian=jacobian)

"""
    successive_relaxation(MVI[, omega=0.])

Solve ``\\vec{x} = \\mathbf{A}^{-1}\\vec{b}`` via the **Successive Relaxation Method**.
Method is Successive Over-Relaxation (SOR) if `omega` > 1, Successive Under-Relaxation (SUR) if `omega` < 1, and reduces to Gauss-Seidel if `omega` = 1.

# Notes
SOR and SUR accelerate or deccelerate convergence of `gauss_seidel()`, respectively, by decreasing or increasing the spectral radius of `MVI.A`.
The core algorithm by which method marches through iterations:
```math
    \\vec{x}^{(k)} = \\bigl( (\\mathbf{D} - ω\\mathbf{L})^{-1} * ((1 - ω)*\\mathbf{D} + ω\\mathbf{U}) \\bigr) ⋅
        \\vec{x}^{(k - 1)} + ω( (\\mathbf{D} - ω\\mathbf{L})^{-1} ) ⋅ \\vec{b}
```

If left unspecified, and if possible, an optimum relaxation parameter, ω will be calculated by `find_omega()`.
"""
successive_relaxation(MVI::MultiVariableIteration, omega::Float64=0.) = solve(MVI; method=:successive_relaxation, omega=omega)

end
