using DataFrames
using LinearAlgebra
using Symbolics

# Ch. 7 (p. 437)
"""
    MultiVariableIteration(f, a, b, n, tol)

Given `f`(p) such that p ∈ [`a`, `b`], find the root of a single-variable equation in so many iterations, `n` within tolerance, `tol`.

---

Iteratively find the solution to a system of equations (SOE): \$\\mathbf{A}\\vec{x} = \\vec{b}\$. Ideal for large, sparse systems.

Parameters
----------
A : tuple
    Either one-dimensional vector of input functions or matrix of characteristic values.
x : tuple
    Either one-dimensional vector of variables or initial guesses for SOE.
b : tuple
    Solution vector.
power : float, optional
    Signed, specified power of tolerance until satisfying method.
max_iter : int, optional
    Number of iterations.
norm_type : {'l_infinity', 'l_two'}, optional
    String representation of desired norm function. `'l_infinity'` by default.

Attributes
----------
A : np.ndarray
    Either one-dimensional vector of input functions or matrix of characteristic values.
x : np.ndarray
    Either one-dimensional vector of variables or initial guesses for SOE.
b : np.ndarray
    Solution vector.
tol : float
    Specified tolerance to which method terminates.
max_iter : int
    Maximum iterations allowed for method.
norm_type : string
    String representation of desired norm function.
is_diagonal, is_symmetric, is_tridiagonal : bool
    Truth value of whether matrix is diagonal, symmetric, and tridiagonal, respectively if not lambda expressions.
eigen_values : np.ndarray
    Eigenvalues of characteristic matrix, `A` if not lambda expressions.
spectral_radius, condition_number : float
    Spectral radius and condition number of characteristic matrix, `A`, respectively if not lambda expressions.

Methods
-------
find_omega(omega=0)
    Suggests optimum \$\\omega\$ over input.
gauss_seidel()
    Improves on `jacobi()` for faster solution.
jacobi()
    Iteratively find solution until within tolerance.
newton_raphson(variables: Tuple[str])
    Given one-dimensional array of equations respect input variables to build gradient (Jacobian) matrix.
successive_relaxation(omega=None)
    Adjusts solution rate of `gauss_seidel()` by scalar \$\\omega\$ which is `None` by default to find the most optimum.

Raises
------
TypeError
    Not all elements in matrix of interest (if one-dimensional) are lambda expressions.
IndexError
    Matrix of interest must be square.
IndexError
    If `x` is not a one-dimensional array.
IndexError
    If `b` is not a one-dimensional array.
ValueError
    If iterations constraint is not an integer greater than zero.
ValueError
    If desired norm method was neither `'l_infinity'` nor `'l_two'`.

See Also
--------
diagonality : Determines if matrix, `A` is strictly, diagonally dominant.
symmetry : Dtermines if matrix, `A` is symmetric.
tridiagonality : Determines if matrix, `A` is tridiagonal.
EigenValues.qr_algorithm : Function to find eigenvalues of matrix, A given initial vector, x and solution vector, b..
spectral_radius : Function to find the spectral radius of characteristic matrix, A.
condition_number : Finds the condition number of matrix, A.
SystemOfEquations : Alternative techniques to solve smaller SOE.

Notes
-----
Specified tolerance evaluated by: `10**power`.

`norm_type` may be either `'l_infinity'` or `'l_two'`. Is 'l_infinity' by default.
"""

struct MultiVariableIteration{T<:Real}
    A   ::Union{Matrix{T}, Vector{Function}}
    x   ::Vector{T}
    b   ::Vector{T}
    N   ::Int64
    tol ::T
end

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
\$\\rho(\\mathbf{A}) = \\max|\\lambda|\$, where \$\\lambda\$ is the set of eigvalsvalues for `A` [burdenNumericalAnalysis2016]_.
"""
spectral_radius(A::Matrix) = maximum(abs.(eigvals(A)))

"""
    condition_number(A)

Find the condition number of matrix, `A`.

# Notes
Definition [burdenNumericalAnalysis2016]_:
    The condition number of the non-singular matrix, \$\\mathbf{A}\$ relative to a norm, \$||\\cdot||\$ is

    \$\$K(\\mathbf{A}) = ||\\mathbf{A}|| \\cdot ||\\mathbf{A}^{-1}||\$\$

A matrix is well-conditioned if \$K(\\mathbf{A})\$ is close to 1 and is ill-conditioned if significantly greater than 1.
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

__find_xk(T, x, c) = T * x + c

"""
    find_omega(mvi[, omega=0.])

Find the optimum relaxation parameter, \$\\omega\$ for `mvi.A` and `mvi.x`, if possible.

# Notes
If 0 < `omega` < 2, then method will converge regardless of choice for `mvi.x`.
If matrix, `mvi.A` is tridiagonal, then the spectral radius of Gauss-Seidel's T-matrix, \$\\mathbf{T}\\_{g}\$ is used to calculate \$\\omega := 2 / (1 + \\sqrt{(1 - `spectral_radius`(\\mathbf{T}\\_{g})})\$ [burdenNumericalAnalysis2016]_.
"""
function find_omega(mvi::MultiVariableIteration, omega::Float64=0.)::Float64
    n, m = size(mvi.A)
    xn, xt = mvi.x, transpose(mvi.x)
    y = (xt * mvi.A) * xn
    theorem_6_22 = issymmetric(mvi.A) && y > 0.
    i, theorem_6_25 = 2, true; while i <= n && theorem_6_25 == true
        theorem_6_25 = det(mvi.A[begin:i, begin:i]) > 0; i += 1
    end
    return if theorem_6_22 || theorem_6_25
        if tridiagonality(mvi.A)
            D   = diagm(diag(mvi.A))
            L   = diagm(-1 => diag(mvi.A, -1))
            U   = diagm(1  => diag(mvi.A, 1))
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

function solve(mvi::MultiVariableIteration;
        method      ::Symbol=:jacobi,
        omega       ::Union{Nothing, Float64}           = nothing,
        variables   ::Union{Nothing, Tuple{Vararg{Num}}}= nothing)
    x = mvi.x
    if method == :newton_raphson
        jacobian = jacobian_form(mvi.A, variables)
        J(x) = convert.(Float64, Symbolics.value.(showjacobian(jacobian, x)))
    elseif method == :successive_relaxation
        if isnothing(omega)
            w = find_omega(mvi)
        # elseif isinstance(omega, (int, float)) && omega > 0.
        #     w = self.find_omega(omega=omega)
        #     logging.info(f"omega = {omega} given. Which is not optimum: {w}")
        #     w = omega
        else
            w = omega
        end
    end
    k, n, approximations, errors = 1, length(x), [x], [mvi.tol*10]
    while errors[end] > mvi.tol && k <= mvi.N
        xk = zeros(eltype(mvi.x), n)
        for i ∈ 1:1:n
            if method == :newton_raphson
                g = [f(x...) for f ∈ mvi.A]
                xk = x + (J(x) \ -g)
            else
                xk[i] = (if method ∈ (:gauss_seidel, :successive_relaxation)
                    y1 = sum([mvi.A[i, j] * xk[j] for j ∈ 1:1:i])
                    y2 = 0.; for j ∈ i+1:1:n
                        y2 += mvi.A[i, j] * x[j]
                    end
                    -y1 - y2
                elseif method == :jacobi
                    -sum([j != i ? mvi.A[i, j] * x[j] : 0. for j ∈ 1:1:n])
                end + mvi.b[i]) / mvi.A[i, i]
            end
        end
        if method == :successive_relaxation
            xk = [(1. - w) * x[i] + w * xk[i] for i ∈ 1:1:n]
        end
        push!(approximations,   xk)
        push!(errors,           norm(xk - x))
        x = xk; k += 1 # iterate to k + 1
    end
    return x
end

"""Given \$\\mathbf{A}\\vec{x} = \\vec{b}\$, use `norm_type` to find \$\\vec{x}\$ via the Gauss-Seidel Method.

Returns
-------
pandas.DataFrame : DataFrame
    Summarized dataframe from iterations.

Attributes
----------
iterations, approximations, errors : np.ndarray
    Collection of iterations, approximations, and normative errors through method.

Warnings
--------
Writes to logfile whether or not a solution was found within the specified tolerance with the supplied, initial guess.

See Also
--------
Norm.l_infinity : Will find \$||x_{i} - x_{0}||_{\\infty}\$
Norm.l_two : Will find \$||x_{i} - x_{0}||_{2}\$

Notes
-----
This improves on `jacobi` by using the most recently calculated entries in the approximation vector, `x` after each iteration.

The primary algorithm by which method marches approximation vector, `x`

.. math::
    \\vec{x}^{(k)} = \\bigl( (\\mathbf{D} - \\mathbf{L})^{-1} * \\mathbf{U} \\bigr) \\cdot \\vec{x}^{(k - 1)} + \\bigl( (\\mathbf{D} - \\mathbf{L})^{-1} \\bigr) \\cdot \\vec{b}
"""
gauss_seidel(mvi::MultiVariableIteration) = solve(mvi; method=:gauss_seidel)

"""Given \$\\mathbf{A}\\vec{x} = \\vec{b}\$, use `norm_type` to find \$\\vec{x}\$ via the Jacobi Method.

Returns
-------
pandas.DataFrame : DataFrame
    Summarized dataframe from iterations.

Attributes
----------
iterations, approximations, errors : np.ndarray
    Collection of iterations, approximations, and normative errors through method.

Warnings
--------
Writes to logfile whether or not a solution was found within the specified tolerance with the supplied, initial guess.

See Also
--------
Norm.l_infinity : Will find \$||x_{i} - x_{0}||_{\\infty}\$
Norm.l_two : Will find \$||x_{i} - x_{0}||_{2}\$

Notes
-----
The primary algorithm by which method marches approximation vector, `x`

.. math::
    \\vec{x}^{(k)} = \\bigl( \\mathbf{D}^{-1} * (\\mathbf{L} + \\mathbf{U}) \\bigr) \\cdot \\vec{x}^{(k - 1)} + ( \\mathbf{D}^{-1} ) \\cdot \\vec{b}
"""
jacobi(mvi::MultiVariableIteration) = solve(mvi; method=:jacobi)

"""Employ the Newton-Raphson Method to find solution of non-linear systems of equations within tolerance.

Parameters
----------
variables : tuple
    Collection of string representations for variables to respect in derivations.

Attributes
----------
iterations, approximations, errors : np.ndarray
    Collection of iterations, approximations, and normative errors through method.

Raises
------
TypeError
    If an element of `variables` is not of type string.

Notes
-----
Modified form of `MultiVariableIteration` to analyze a one-dimensional array of non-linear SOE. Each element should be a lambda expression wherein each variable is represented.

Examples 
--------
>>> A = [lambda x1, x2, x3: 3*x1 - sympy.cos(x2*x3) - 1/2,
    lambda x1, x2, x3: x1**2 - 81*(x2 + 0.1)**2
        + sympy.sin(x3) + 1.06,
    lambda x1, x2, x3: sympy.exp(-x1*x2)
        + 20*x3 + (10*math.pi - 3)/3
    ]
>>> x, b = (0.1, 0.1, -0.1), (0, 0, 0)
>>> variables = ("x1", "x2", "x3")
>>> MultiVariableIteration(A, x, b).newton_raphson(variables)["Approximations"].values[-1]
[0.5, 0., -0.52359877]
"""
newton_raphson(mvi::MultiVariableIteration, variables::Tuple{Vararg{Num}}) = solve(mvi;
method=:newton_raphson, variables=variables)

"""Given \$\\mathbf{A}\vec{x} = \\vec{b}\$, use `norm_type` to find \$\\vec{x}\$ via the Successive Relaxation Method. Is Successive Over-Relaxation (SOR) if `omega` > 1, Successive Under-Relaxation (SUR) if `omega` < 1, and is Gauss-Seidel if `omega` = 1.

Parameters
----------
omega : None or float, optional
    Relaxation parameter.

Attributes
----------
iterations, approximations, errors : np.ndarray
    Collection of iterations, approximations, and normative errors through method.

Returns
-------
pandas.DataFrame : DataFrame
    Summarized dataframe from iterations.

Warnings
--------
Writes to logfile optimal choice of omega, regardless of assignment, and whether or not a solution was found within the specified tolerance with the supplied, initial guess.

See Also
--------
find_omega : Will analyze SOE to find an optimal \$\\omega\$, if possible.
gauss_seidel : Gauss-Seidel Method modified by omega.
Norm.l_infinity : Will find \$||x_{i} - x_{0}||_{\\infty}\$
Norm.l_two : Will find \$||x_{i} - x_{0}||_{2}\$

Notes
-----
SOR and SUR modify, respectively, on `gauss_seidel` by decreasing or increasing, respectively, the spectral radius of `A` to accelerate or deccelerate convergence, respectively.

The primary algorithm by which method marches approximation vector, `x`

.. math::
    \\vec{x}^{(k)} = \\bigl( (\\mathbf{D} - \\omega\\mathbf{L})^{-1} * ((1 - \\omega)*\\mathbf{D} + \\omega\\mathbf{U}) \\bigr) \\cdot \\vec{x}^{(k - 1)} + \\omega( (\\mathbf{D} - \\omega\\mathbf{L})^{-1} ) \\cdot \\vec{b}

which is similar to `gauss_seidel`

.. math::
    \\vec{x}^{(k)} = \\bigl( (\\mathbf{D} - \\mathbf{L})^{-1} * \\mathbf{U} \\bigr) \\cdot \\vec{x}^{(k - 1)} + \\bigl( (\\mathbf{D} - \\mathbf{L})^{-1} \\bigr) \\cdot \\vec{b}

`omega` will be analyzed independent of assigned value which will be used if not specified in assignment and if possible.
"""
successive_relaxation(mvi::MultiVariableIteration, omega::Union{Nothing, Float64}=nothing) = solve(mvi; method=:successive_relaxation, omega=omega)