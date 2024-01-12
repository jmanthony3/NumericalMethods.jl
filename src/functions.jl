using IntervalArithmetic, IntervalRootFinding
using Polynomials
using Printf
using StaticArrays
using Statistics
using Symbolics

# Ch. 2 (p. 47)
"""
    SVI(f, a, b, n[, tol=10^-3; g])

Given `f`(`x`) such that `x` ∈ [`a`, `b`], find the root of a single-variable, equation in so many iterations, `n` within tolerance, `tol`.

Methods
-------
bisection()
    Search for solution by halving the bounds wherein `a` and `b` initially yield opposite signs in function.
false_position(p0: float, p1: float)
    solution bounded by `a` and `b` wherein initial guesses `p0` and `p1` yield opposite signs in function.
fixed_point(p0: float)
    Root-finding method to find solution near initial guess.
newton_raphson(p0: float)
    Root-finding method to find solution near initial guess.
secant_method(p0: float, p1: float)
    Initial guesses `p0` and `p1` must yield opposite signs in function. Solution is NOT bounded by `a` and `b`.

# Notes
Convergence Rates:

    `newton_raphson` > `secant_method` > `false_position` > `fixed_point` > `bisection`
"""
struct SVI
    f   ::Function
    a   ::Real
    b   ::Real
    n   ::Integer
    tol ::Real
end

"""
    maximumslope(obj::SVI)

Find the greatest value for first derivative of function.
"""
function maximumslope(obj::SVI)::AbstractFloat
    @variables x
    Dx  = Differential(x)
    df  = simplify(expand_derivatives(Dx(obj.f(x))); expand=true)
    df  = build_function(df, x, expression=Val{false})
    return maximum(abs.(df.(range(obj.a, obj.b, length=1000))))
end

"""
    max_iterations(obj, method[, p0=0])

Find greatest integer for maximum iterations within tolerance.

# Notes
Acceptable values for `method` ∈ {`:bisection`, `:fixed_point`, `:newton_raphson`, `:secant_method`, `:false_position`}.
Initial guess, `p0` for function solution is not needed for `:bisection` method.
"""
function max_iterations(
    obj     ::SVI,
    method  ::Symbol,
    p0      ::Real      = 0.;
    k       ::Real      = NaN
)::Integer
    if method == "bisection"
        return ceil(-log(obj.tol / (obj.b - obj.a))
            / log(2))
    elseif method in ("fixed_point", "newton_raphson", "secant_method", "false_position")
        return ceil(-log(obj.tol / max(p0 - obj.a, obj.b - p0))
            / log(k))
    else
        error("Method must be: `:bisection`, `:fixed_point`, `:newton_raphson`, `:secant_method`, or `:false_position`.")
    end
    # logging.info(f"With the inputs, I will terminate the technique after so many iterations, N = {max_iter}")
end

## 2.1 (p. 48)
"""
    bisection(obj::SVI)

Root-finding method: f(x) = 0.
Search for solution by halving the bounds wherein `a` and `b` initially yield opposite signs in function.

# Notes
Relying on the Intermediate Value Theorem (IVT), this is a bracketed, root-finding method.
This method is rather slow to converge but will always converge to a solution; therefore, is a good starter method.
"""
function bisection(obj)
    f, a, b = obj.f, obj.a, obj.b
    if f(a)*f(b) < 0    # check if f(a) and f(b) are opposite signs
        N = obj.n
        # initialize
        k = 1
        g = zeros(MVector{N})
        r = zeros(MVector{N})
        g[k], r[k] = f(a), 1.
        # exit by whichever condition is `true` first
        while r[k] >= obj.tol && k < N
            x = (b - a) / 2.
            p = a + x                           # new value, p
            f(a)*f(p) > 0 ? (a = p) : (b = p)   # adjust next bounds
            g[k+1], r[k+1] = p, abs(x)          # error of new value, p
            k += 1                              # iterate to k + 1
        end
        return k <= N ? g[k] : NaN
    else                # abort if f(a) is not opposite f(b)
        error(@sprintf("Interval bounds must yield opposite signs in function, f := [f(a = %1.4f) = %1.4f, f(b = %1.4f) = %1.4f]",
            a, b, f(a), f(b)))
    end
end

## 2.2 (p. 55)
"""
    fixed_point(obj, p0)

Attempt root-finding method with initial guess, `p0` in [a, b] by solving the equation g(p) = p via f(p) - p = 0.

**Use function with lowest slope!**

_Not root-bracketed._

# Notes
Theorem:
1) Existence of a fixed-point:
    If g ∈ C[a,b] and g(x) ∈ C[a, b] for all x ∈ [a, b], then function, g has a fixed point, p ∈ [a, b].
2) Uniqueness of a fixed point:
    If g'(x) exists on [a, b] and a positive constant, k < 1 exist with {|g'(x)| ≤ k | x ∈ (a, b)}, then there is exactly one fixed-point, p ∈ [a, b].

Converges by mathcal{O}(text{linear}) if g'(p) ≠ 0, and mathcal{O}(text{quadratic}) if g'(p) = 0 and g''(p) < M, where M = g''(ξ) that is the error function.
"""
function fixed_point(obj::SVI, p0::AbstractFloat)::AbstractFloat
    f, a, b = obj.f, obj.a, obj.b
    N = obj.n
    # initialize
    k = 1
    g = zeros(MVector{N})
    r = zeros(MVector{N})
    g[k], r[k] = f((a + b) / 2.), 1.
    # exit by whichever condition is `true` first
    while r[k] >= obj.tol && k < N
        p = f(p0)                           # new value, p
        g[k+1], r[k+1] = p, abs((p-p0)/p0)  # error of new value, p
        p0 = p; k += 1                      # iterate to k + 1
    end
    return k <= N ? g[k] : NaN
end

## 2.3 (p. 66)
"""
    newton_raphson(obj, p0)

Attempt root-finding method with initial guess, `p0` in [a, b] by solving the equation g(p) = p via f(p) - p = 0.

**Use function with lowest slope!**

_f'(x) ≠ 0_

# Notes
Not root-bracketed and has trouble with symmetric functions!
Initial guess, `p0` must be close to real solution; else, will converge to different root or oscillate (if symmetric).
This method can be viewed as fixed-point iteration.

Technique based on first Taylor polynomial expansion of f about p_{0} (that is `p0`) and evaluated at x = p. |p - p_{0}| is assumed small; therefore, 2^{\text{nd}}-order Taylor term, the error, is small.

Newton-Raphson has quickest convergence rate.

See `fixed_point()` for theorem.
"""
function newton_raphson(obj::SVI, p0::AbstractFloat)::AbstractFloat
    f, a, b = obj.f, obj.a, obj.b
    # determine form of derivative
    @variables x
    Dx  = Differential(x)
    df  = simplify(expand_derivatives(Dx(f(x))); expand=true)
    df  = build_function(df, x, expression=Val{false})
    N = obj.n
    # initialize
    k = 1
    g = zeros(MVector{N})
    r = zeros(MVector{N})
    g[k], r[k] = f(a), 1.
    # exit by whichever condition is TRUE first
    while r[k] >= obj.tol && k < N
        p = p0 - (f(p0) / df(p0))           # new value, p
        g[k+1], r[k+1] = p, abs(p - p0)     # error of new value, p
        p0 = p; k += 1                      # iterate to k + 1
    end
    return k <= N ? g[k] : NaN
end

"""
    secant_method(obj, p0, p1)

Attempt root-finding method with initial guesses, `p0` and `p1` in [a, b] by solving the equation g(p) = p via f(p) - p = 0.

**Use function with lowest slope!**

_Not root-bracketed._

# Notes
Method is less computationally expensive than `newton_raphson()` but may converge at slower rate by circumventing need to calculate derivative.

See `fixed_point()` for theorem.
"""
function secant_method(obj::SVI, p0::T, p1::T)::AbstractFloat where {T<:AbstractFloat}
    f, a, b = obj.f, obj.a, obj.b
    if f(p0)*f(p1) < 0  # check if f(p0) and f(p1) are opposite signs
        N = obj.n
        # initialize
        k = 1
        g = zeros(MVector{N})
        r = zeros(MVector{N})
        g[k], r[k] = f(a), 1.
        # exit by whichever condition is TRUE first
        while r[k] >= obj.tol && k < N
            q0, q1 = f(p0), f(p1)
            p = p1 - q1*(p1 - p0)/(q1 - q0)     # new value, p
            g[k+1], r[k+1] = p, (abs(p - p0))   # error of new value
            p0, p1 = p1, p; k += 1              # iterate to k + 1
        end
        return k <= N ? g[k] : NaN
    else                # abort if f(p0) is not opposite f(p1)
        error(@sprintf("Interval bounds must yield opposite signs in function, f := [f(p0 = %1.4f) = %1.4f, f(p1 = %1.4f) = %1.4f]",
            p0, p1, f(p0), f(p1)))
    end
end

"""
    false_position(obj, p0, p1)

Attempt root-finding method with initial guesses, `p0` and `p1` in [a, b] solving the equation g(p) = p via f(p) - p = 0.

**Use function with lowest slope!**

# Notes
Similar to, but slower to converge than, the `secant_method()` by including a test to ensure solution is root-bracketed.

See `fixed_point()` for theorem.
"""
function false_position(obj::SVI, p0::T, p1::T)::AbstractFloat where {T<: AbstractFloat}
    f, a, b = obj.f, obj.a, obj.b
    if f(p0)*f(p1) < 0  # check if f(p0) and f(p1) are opposites signs
        N = obj.n
        # initialize
        k = 1
        g = zeros(MVector{N})
        r = zeros(MVector{N})
        g[k], r[k] = f(a), 1.
        # exit by whichever condition is `true` first
        while r[k] >= obj.tol && k < N
            q0, q1 = f(p0), f(p1)
            p = p1 - q1*(p1 - p0)/(q1 - q0)     # new value, p
            g[k+1], r[k+1] = p, abs(p - p0)     # error of new value, p
            f(p)*q1 < 0 ? p0 = p1 : nothing     # adjust next bounds
            p1 = p; k += 1                      # iterate to k + 1
        end
        return k <= N ? g[k] : NaN
    else                # abort if f(p0) is not opposite f(p1)
        error(@sprintf("Interval bounds must yield opposite signs in function, f := [f(p0 = %1.4f) = %1.4f, f(p1 = %1.4f) = %1.4f]",
            p0, p1, f(p0), f(p1)))
    end
end

# Ch. 3 (p. 103)
## 3.1 (p. 104)
"""
    lagrange(x, f[; n=nothing])

Given a domain, `x` and range, `f`, construct the `n`th Lagrangian polynomial.

# Notes
If `n=nothing`, then method will utilize entire dataset.
Polynomials will quickly oscillate for larger datasets.
"""
function lagrange(
    x::T,
    f::T;
    n::Union{Integer, Nothing}   = nothing
)::Tuple{Function, AbstractVector} where {T<:AbstractVector}
    @variables t
    Dt = Differential(t)
    function coefficient(xₖ, t)
        num, den = [], []
        for xₗ ∈ x
            if isa(xₗ, Num) || xₗ != xₖ
                push!(num, (t - xₗ))
                push!(den, (xₖ - xₗ))
            end
        end
        return prod(num) / prod(den)
    end
    function error(n, ξ, t)
        s       = []
        ξ_error = zeros(MVector{n})
        for i ∈ 1:1:n
            push!(s, t - x[i])
            g               = simplify(prod(s); expand=true)
            g_eval          = build_function(g, t, expression=Val{false})
            gp              = simplify(expand_derivatives(Dt(g)); expand=true)
            gp_eval         = build_function(gp, t, expression=Val{false})
            coeffs_sorted   = if hasproperty(Symbolics.value(gp), :dict)
                coeffs_map      = Symbolics.value(gp).dict
                exps            = Symbolics.degree.(keys(coeffs_map))
                collect(values(coeffs_map))[sortperm(exps)]
            else
                []
            end
            if gp_eval(0.) != 0.
                pushfirst!(coeffs_sorted, gp_eval(0.))
            end
            g_prime_poly    = Polynomial(float.(coeffs_sorted))
            gx              = if i == 1
                g_eval(coeffs_sorted[begin])
            elseif i == 2
                g_eval(roots(g_prime_poly)[begin])
            else
                R = []
                for r ∈ roots(g_prime_poly)
                    if isreal(r)
                        push!(R, g_eval(r))
                    end
                end
                maximum(abs.(R))
            end
            for k ∈ 1:1:i
                ξ = simplify(expand_derivatives(Dt(ξ)); expand=true)
            end
            dξ              = maximum(
                build_function(ξ, t; expression=Val{false}).(x[begin:i])
                    ./ (factorial(n)))
            ξ_error[i]      = dξ * abs(gx)
        end
        return maximum(abs.(ξ_error))
    end
    n  = (isnothing(n) ? length(x) - 1 : n)
    terms   = []
    errors  = zeros(MVector{n + 1})
    for k ∈ 1:1:n + 1
        push!(terms, f[k] * coefficient(x[k], t))
        errors[k]   = error(k, sum(terms), t)
    end
    p       = build_function(sum(terms), t, expression=Val{false})
    return p, errors
end

## 3.3 (p. 122)
"""
    newtondifference(x, f, α[; dir::Symbol=:auto])

Given a domain, `x` and range, `f`, construct some polynomial by Newton's Divided Difference centered around `α`.
`:forward` or `:backward` construction.

# Notes
Direction will be chosen if not specified.
Polynomials best made with even spacing in `x`; although, this is not completely necessary.
"""
function newtondifference(
    x   ::T,
    f   ::T,
    α   ::Real;
    dir ::Symbol        = :auto
)::Function where {T<:AbstractVector}
    if dir == :auto
        dir = (α <= median(x) ? :forward : :backward)
    end
    fterm(g, i ,j) = (g[i, j] - g[i - 1, j]) / (g[i, 1] - g[i - (j - 1), 1])
    m, n    = length(x), length(x) + 1
    fxn     = MMatrix{m, n}(zeros((m, n)))
    coeff   = zeros(MVector{m - 1})
    fxn[:, 1], fxn[:, 2] = x, f
    for j ∈ 2:1:m
        for i ∈ j:1:m
            fₖ = fterm(fxn, i, j)
            fxn[i, j + 1] = fₖ
            if dir == :forward && i == j
                coeff[j - 1] = fₖ
            elseif dir == :backward && i == m
                coeff[j - 1] = fₖ
            end
        end
    end
    @variables t
    k, g, terms = (dir == :forward ? 1 : m), 0., 1.
    for c ∈ coeff
        terms  *= (t - x[k])
        g      += c * prod(terms)
        k      += (dir == :forward ? 1 : -1)
    end
    p = g + (dir == :forward ? f[begin] : f[end])
    return build_function(p, t, expression=Val{false})
end

## 3.5 (p. 142)
"""
    natural(x, f)

The bookend polynomials do not assume the slope entering and exiting the interval as the derivative at the respective endpoint.
"""
function natural(
    x   ::T,
    f   ::T
)::Tuple{AbstractVector, AbstractArray{Function}} where {T<:AbstractVector}
    function _algorithm(g)
        Y               = g
        m, n            = length(Y), length(Y) - 1
        # STEP 1:   build list, h_i
        H               = zeros(MVector{n})
        for i ∈ 1:1:n
            H[i] = X[i+1] - X[i]
        end
        # STEP 2:   build list, alpha_i
        A, ALPHA        = Y, zeros(MVector{m})
        # ALPHA[1]        = 3*(A[2] - A[1])/H[1] - 3*AP[1]
        # ALPHA[m]        = 3*AP[m] - 3*(A[m] - A[n])/H[n]
        for i ∈ 2:1:n
            ALPHA[i] = 3/H[i]*(A[i+1] - A[i]) - 3/H[i-1]*(A[i] - A[i-1])
        end
        # Algorithm 6.7 to solve tridiagonal
        # STEP 3:   define l, mu, and z first points
        L, MU           = zeros(MVector{m}), zeros(MVector{m})
        Z, C            = zeros(MVector{m}), zeros(MVector{m})
        L[1], MU[1], Z[1]= 1., 0, 0.
        # STEP 4:   build lists l, mu, and z
        for i ∈ 2:1:n
            L[i]  = 2(X[i+1] - X[i-1]) - H[i-1]*MU[i-1]
            MU[i] = H[i] / L[i]
            Z[i]  = (ALPHA[i] - H[i-1]*Z[i-1]) / L[i]
        end
        # STEP 5:   define l, z, and c endpoints
        L[m], Z[m], C[m]= 1., 0., 0.
        # STEP 6:   build lists c, b, and d
        B, D            = zeros(MVector{n}), zeros(MVector{n})
        for i ∈ 0:1:n-1
            j    = n-i
            C[j] = Z[j] - MU[j]*C[j+1]
            B[j] = (A[j+1] - A[j])/H[j] - H[j]*(C[j+1] + 2C[j])/3
            D[j] = (C[j+1] - C[j]) / 3H[j]
        end
        return Y, A, B, C, D
    end
    g, X            = f, x
    Y, A, B, C, D   = _algorithm(g)
    n, splines      = length(X) - 1, []
    for j ∈ 1:1:n
        xj, aj, bj, cj, dj = X[j], A[j], B[j], C[j], D[j]
        sj(x) = aj + bj*(x - xj) + cj*(x - xj)^2 + dj*(x - xj)^3
        push!(splines, sj)
    end
    return Y, splines
end

"""
    clamped(x, f, fp)

The bookend polynomials will have the same slope entering and exiting the interval as the derivative at the respective endpoint.
"""
function clamped(
    x   ::T,
    f   ::T,
    fp  ::T
)::Tuple{AbstractVector, AbstractArray{Function}} where {T<:AbstractVector}
    function _algorithm(g, gp)
        Y, YP           = g, gp
        m, n            = length(Y), length(Y) - 1
        # STEP 1:   build list, h_i
        H               = zeros(MVector{n})
        for i ∈ 1:1:n
            H[i] = X[i+1] - X[i]
        end
        # STEP 2:   define alpha list endpoints
        A, AP, ALPHA    = Y, YP, zeros(MVector{m})
        ALPHA[1]        = 3(A[2] - A[1])/H[1] - 3AP[1]
        ALPHA[m]        = 3AP[m] - 3(A[m] - A[n])/H[n]
        # STEP 3:   build list, alpha_i
        for i ∈ 2:1:n
            ALPHA[i] = 3/H[i]*(A[i+1] - A[i]) - 3/H[i-1]*(A[i] - A[i-1])
        end
        # Algorithm 6.7 to solve tridiagonal
        # STEP 4:   define l, mu, and z first points
        L, MU           = zeros(MVector{m}), zeros(MVector{m})
        Z, C            = zeros(MVector{m}), zeros(MVector{m})
        L[1], MU[1]     = 2H[1], 0.5
        Z[1]            = ALPHA[1] / L[1]
        # STEP 5:   build lists l, mu, and z
        for i ∈ 2:1:n
            L[i]  = 2(X[i+1] - X[i-1]) - H[i-1]*MU[i-1]
            MU[i] = H[i]/L[i]
            Z[i]  = (ALPHA[i] - H[i-1]*Z[i-1])/L[i]
        end
        # STEP 6:   define l, z, and c endpoints
        L[m]            = H[n] * (2 - MU[n])
        Z[m]            = (ALPHA[m] - H[n]*Z[n]) / L[m]
        C[m]            = Z[m]
        # STEP 7:   build lists c, b, and d
        B, D            = zeros(MVector{n}), zeros(MVector{n})
        for i ∈ 0:1:n-1
            j    = n-i
            C[j] = Z[j] - MU[j]*C[j+1]
            B[j] = (A[j+1] - A[j])/H[j] - H[j]*(C[j+1] + 2*C[j])/3
            D[j] = (C[j+1] - C[j]) / 3H[j]
        end
        return Y, A, B, C, D
    end
    g, X, gp        = f, x, fp
    Y, A, B, C, D   = _algorithm(g, gp)
    n, splines      = length(X) - 1, []
    for j ∈ 1:1:n
        xj, aj, bj, cj, dj = X[j], A[j], B[j], C[j], D[j]
        sj(x) = aj + bj*(x - xj) + cj*(x - xj)^2 + dj*(x - xj)^3
        push!(splines, sj)
    end
    return Y, splines
end

## 3.6 (p. 162)
"""
    bezier(x, y, xguides, yguides)

An application of Hermitic polynomials to draw Bezier curves between points.

# Notes
Each argument should be a one-to-one mapping of points, (xᵢ, yᵢ) and (xᵢ₊₁, yᵢ₊₁) and their respective guide points, (xᵢ⁺, yᵢ⁺) and (xᵢ₊₁⁻, yᵢ₊₁⁻).
"""
function bezier(x::T, y::T, xguides::T, yguides::T)::AbstractArray{Function} where {T<:AbstractVector}
    n, curves = length(x) - 1, []
    for i ∈ 1:1:n
        a = (x[i],
            3(xguides[i] - x[i]),
            3(x[i] + xguides[i + 1] - 2xguides[i]),
            x[i + 1] - x[i] + 3xguides[i] - 3xguides[i + 1])
        b = (y[i],
            3(yguides[i] - y[i]),
            3(y[i] + yguides[i + 1] - 2yguides[i]),
            y[i + 1] - y[i] + 3yguides[i] - 3yguides[i + 1])
        xi(t) = sum(a .* (1., t, t^2, t^3))
        yi(t) = sum(b .* (1., t, t^2, t^3))
        push!(curves, (x = xi, y = yi))
    end
    return curves
end

# Ch. 4 (p. 171)
## 4.1 (p. 172)
"""
    n1derivative(x, f, j[; n=nothing])

The general (n + 1)-point formula to approximate f' at point `j`.

# Notes
If `n = nothing`, then entire dataset used to construct `n`th Lagrange coefficient.
"""
function n1derivative(
    x::T,
    f::T,
    j::Integer;
    n::Union{Integer, Nothing}   = nothing
)::AbstractFloat where {T<:AbstractVector}
    @variables t
    Dt = Differential(t)
    function coefficient(xₖ, t)
        num, den = [], []
        for xₗ ∈ x
            if isa(xₗ, Num) || xₗ != xₖ
                push!(num, (t - xₗ))
                push!(den, (xₖ - xₗ))
            end
        end
        return prod(num) / prod(den)
    end
    n  = (isnothing(n) ? length(x) - 1 : n)
    gp      = 0.
    for k ∈ 1:1:n + 1
        Lₖ          = coefficient(x[k], t)
        Lₖp         = simplify(expand_derivatives(Dt(Lₖ)); expand=true)
        Lₖp_eval    = build_function(Lₖp, t, expression=Val{false})
        gp         += f[k]*Lₖp_eval(x[j])
    end
    return gp
end

"""
    endpoint(x, f, h, point[; method=:three])

Find the derivative of a bookend point at either `:begin` or `:end` of dataset.
Acceptable values for `method` include {`:three`, `:five`}.
"""
function endpoint(
    x       ::T,
    f       ::T,
    h       ::Real,
    point   ::Symbol;
    method  ::Symbol    = :three
)::AbstractFloat where {T<:AbstractVector}
    i = (point == :begin ? 1 : (point == :end ? length(x) : nothing))
    if point == :begin
        if method == :three
            (-3f[i] + 4f[i+1] - f[i+2]) / 2h
        elseif method == :five
            (-25f[i] + 48f[i+1]
                - 36f[i+2] + 16f[i+3]
                    - 3f[i+4]) / 12h
        else
            NaN
        end
    elseif point == :end
        if method == :three
            -(-3f[i] + 4f[i-1] - f[i-2]) / 2h
        elseif method == :five
            -(-25f[i] + 48f[i-1]
                - 36f[i-2] + 16f[i-3]
                    - 3f[i-4]) / 12h
        else
            NaN
        end
    else
        NaN
    end
end

"""
    midpoint(x, f, h, point[; method=:three])

Find the derivative of some point within a dataset.
Acceptable values for `method` include {`:three`, `:five`, `:2nd`}.
"""
function midpoint(
    x       ::T,
    f       ::T,
    h       ::Real,
    point   ::Integer;
    method  ::Symbol    = :three
)::AbstractFloat where {T<:AbstractVector}
    i = point
    return if method == :three
        (f[i + 1] - f[i - 1]) / 2h
    elseif method == :five
        (f[i - 2] - 8f[i - 1]
            + 8f[i + 1] - f[i + 2]) / 12h
    elseif method == :2nd
        (f[i - 1] - 2f[i] + f[i + 1]) / (h ^ 2.)
    else
        NaN
    end
end

## 4.3 (p. 191)
"""
    integrate(f, x      [; rule=:trapezoidal, tol=10^-3])
    integrate(f, a, b, h[; rule=:trapezoidal, tol=10^-3])
    integrate(f, a, b, n[; rule=:trapezoidal, tol=10^-3])

Find the definite integral by some numerical quadrature.

# Notes
`f` may be a function or range.
The domain may be defined with a vector, `x` or on the interval [`a`, `b`] either by number of sub-intervals, `n` or step-size, `h`.
`rule` accepts {`:trapezoidal` (default), `:midpoint`, `:simpson13`, `:simpson38`, `:simpsonN`}.
Dataset may contain unevenly spaces points.

# References
https://en.wikipedia.org/wiki/Simpson%27s_rule
"""
function integrate(
    f       ::Union{AbstractVector, Function},
    x       ::AbstractVector;
    rule    ::Symbol    = :trapezoidal,
    tol     ::Real      = 10^-3
)::AbstractFloat
    is_function = isa(f, Function)
    a, b, n     = x[begin], x[end], length(x) - 1
    if is_function
        n_min = if rule == :trapezoidal
            ceil(sqrt((b - a)^3 / (12 * tol)))
        elseif rule == :simpson13
            ceil(((b - a)^5 / (180 * tol)) ^ (1/4))
        elseif rule == :midpoint
            ceil(sqrt((b - a)^3 / (6 * tol)))
        end
    end
    if rule ∈ [:simpson13, :midpoint] && isodd(n)
        F           = integrate(is_function ? f : f[n:end], x[n:end], rule=:trapezoidal)
        x           = x[begin:n]
        a, b, n     = x[begin], x[end], length(x) - 1
    elseif rule == :simpson38 && n % 3 != 0
        m           = n - (n % 3 - 1)
        F           = integrate(is_function ? f : f[m:end], x[m:end], rule=:trapezoidal)
        x           = x[begin:m]
        a, b, n     = x[begin], x[end], length(x) - 1
    elseif rule == :simpsonN && isodd(n)
        hn2, hn1    = x[n] - x[n - 1], x[n + 1] - x[n]
        fn2, fn1, fn = is_function ? f.(x[n - 1:n + 1]) : f[n - 1:n + 1]
        alpha       = (2hn1^2 + 3hn1*hn2) / 6(hn2 + hn1)
        beta        = (hn1^2 + 3hn1*hn2) / 6hn2
        eta         = (hn1^3) / (6hn2 * (hn2 + hn1))
        F           = alpha*fn + beta*fn1 - eta*fn2
    else
        F           = 0.
    end
    if rule == :trapezoidal
        h           = (b - a) / n
        z           = 0.
        for j ∈ 2:1:n
            z += (is_function ? f(x[j]) : f[j])
        end
        F          += if is_function
            h/2*(f(a) + 2z + f(b))
        else
            h/2*(f[begin] + 2z + f[end])
        end
    elseif rule == :simpson13
        h           = (b - a) / n
        z1          = 0.
        for j ∈ 2:1:(n ÷ 2)
            z1 += (is_function ? f(x[2j - 1]) : f[2j - 1])
        end
        z2          = 0.
        for j ∈ 1:1:(n ÷ 2)
            z2 += is_function ? f(x[2j]) : f[2j]
        end
        F          += if is_function
            h/3*(f(a) + 2z1 + 4z2 + f(b))
        else
            h/3*(f[begin] + 2z1 + 4z2 + f[end])
        end
    elseif rule == :simpson38
        h           = (b - a) / n
        z1          = 0.
        for j ∈ 2:1:n
            if j % 3 != 0
                z1 += (is_function ? f(x[j]) : f[j])
            end
        end
        z3          = 0.
        for j ∈ 1:1:(n ÷ 3)
            z3 += (is_function ? f(x[3j]) : f[3j])
        end
        F          += if is_function
            3h/8*(f(a) + 3z1 + 2z3 + f(b))
        else
            3h/8*(f[begin] + 3z1 + 2z3 + f[end])
        end
    elseif rule == :simpsonN
        h           = (b - a) / n
        for j ∈ 0:1:(n ÷ 2) - 1
            h2j0, h2j1 = x[2j + 2] - x[2j + 1], x[2j + 3] - x[2j + 2]
            f2j0, f2j1, f2j2 = is_function ? f.(x[2j + 1:2j + 3]) : f[2j + 1:2j + 3]
            F += 6 \ (h2j0 + h2j1) * (
                (2 - h2j1 / h2j0) * f2j0
                    + ((h2j0 * h2j1) \ (h2j0 + h2j1)^2) * f2j1
                    + (2 - h2j0 / h2j1) * f2j2)
        end
    elseif rule == :midpoint
        h           = (b - a) / (n + 2)
        z           = 0.
        for j ∈ 1:1:(n ÷ 2)
            z += (is_function ? f(x[2j]) : f[2j])
        end
        F          += 2h*z
    end
    return F
end

function integrate(
    f       ::Union{AbstractVector, Function},
    a       ::Real,
    b       ::Real,
    h       ::Real;
    rule    ::Symbol    = :trapezoidal,
    tol     ::Real      = 10^-3
)::AbstractFloat
    return integrate(f, float.(a:h:b), rule=rule, tol=tol)
end

function integrate(
    f       ::Union{AbstractVector, Function},
    a       ::Real,
    b       ::Real,
    n       ::Integer;
    rule    ::Symbol    = :trapezoidal,
    tol     ::Real      = 10^-3
)::AbstractFloat
    return if rule == :midpoint
        integrate(f, float.(a:(b - a)/(n + 2):b), rule=rule, tol=tol)
    else
        integrate(f, float.(a:(b - a)/n:b), rule=rule, tol=tol)
    end
end

# Ch. 5 (p. 259)
"""
    ODE(f, a, b, h, α, β, N)

Structure of the boundary conditions to differential equation.

# Notes
Make sure the independent variable (e.g. time) is the first argument of `f`!
"""
struct ODE
    f::Function
    a::Real
    b::Real
    h::Real
    α::Real
    β::Real
    N::Integer
end

"""
    ivp(obj::ODE[; tol=10^-3, method=:forward_euler])

Solve `obj` according to `method` ∈ {`:forward_euler` (default), `:backward_euler`, `:improved_euler`, `:modified_euler`, `:runge_kutta`}.

# Notes
Each method has an equivalent convenience function.
E.g. `ivp(obj; method=:runge_kutta)` ≡ `runge_kutta(obj)`.
"""
function ivp(obj::ODE; tol=10^-3, method=:forward_euler)::AbstractVector
    f       = obj.f
    t, h, w = obj.a, obj.h, obj.α
    ea, eb, λ = 1/2, 1/2, 1
    g       = zeros(MVector{obj.N + 1})
    g[1]    = w
    for i ∈ 1:1:obj.N
        w += if method == :forward_euler
            h * f(t, w)
        elseif method == :backward_euler
            h * f(t + h, w + h*f(t, w))
        elseif method ∈ [:improved_euler, :modified_euler]
            h * (ea*f(t, w) + eb*f(t + λ*h, w + λ*h*f(t, w)))
        elseif method == :runge_kutta
            k1 = h * f(t,       w)
            k2 = h * f(t + h/2, w + k1/2)
            k3 = h * f(t + h/2, w + k2/2)
            k4 = h * f(t + h,   w + k3)
            (k1 + 2k2 + 2k3 + k4) / 6
        end
        g[i + 1] = w
        # push!(increments, w - w0)
        t   = obj.a + i*h
    end
    return g
end

forward_euler(obj;  tol=10^-3) = ivp(obj; tol=tol, method=:forward_euler)
backward_euler(obj; tol=10^-3) = ivp(obj; tol=tol, method=:backward_euler)
improved_euler(obj; tol=10^-3) = ivp(obj; tol=tol, method=:improved_euler)
modified_euler(obj; tol=10^-3) = ivp(obj; tol=tol, method=:modified_euler)
runge_kutta(obj;    tol=10^-3) = ivp(obj; tol=tol, method=:runge_kutta)

# Ch. 8 (p. 505)
## 8.1 (p. 506)
"""
    linearleastsquares(x, f, n::Integer)

Construct a polynomial of degree, `n` while minimizing the least squares error.

# Notes
Least squares error := ``E = \\sum_{i=1}^{m}[y_{i} - P_{n}(x_{i})]^{2}``

Constructed polynomial of the form: ``P(x) = a_{n}x^{n} + a_{n - 1}x^{n - 1} + \\dots + a_{1}x + a_{0}``
"""
function linearleastsquares(
    x::T,
    f::T,
    n::Integer
)::Tuple{Function, AbstractFloat} where {T<:AbstractVector}
    X, Y    = x, f
    m       = length(X)
    A       = MMatrix{n+1, n+1}(zeros((n+1, n+1)))
    b       = zeros(MVector{n+1})
    for i ∈ 0:1:n
        for j ∈ 0:1:n
            for k ∈ 1:1:m
                A[i+1,j+1] += X[k]^(i + j)
            end
        end
        for j ∈ 1:1:m
            b[i+1] += Y[j]*X[j]^i
        end
    end
    @variables t
    x, p    = A\b, 0.
    for i ∈ 0:1:n
        p  += x[i+1]*t^i
    end
    p       = build_function(p, t, expression=Val{false})
    error   = sum((Y - p.(X)) .^ 2)
    return p, error
end

"""
    linearleastsquares(x, f, type::Symbol)

Given a domain and range, yield the coefficients for an equation and the equation of the form ``y = ax^{b}``.
"""
function linearleastsquares(
    x       ::T,
    f       ::T,
    type    ::Symbol
)::Tuple{Function, AbstractFloat} where {T<:AbstractVector}
    if type == :power
        X, Y    = x, f
        m       = length(X)
        q1      = zeros(MVector{m})
        q2      = zeros(MVector{m})
        q3      = zeros(MVector{m})
        q4      = zeros(MVector{m})
        for i ∈ 1:1:m
            q1[i] = log(X[i])*log(Y[i])
            q2[i] = log(X[i])
            q3[i] = log(Y[i])
            q4[i] = log(X[i])^2
        end
        num     = m*sum(q1) - sum(q2)*sum(q3)
        den     = m*sum(q4) - sum(q2)^2
        b       = num / den
        a       = exp((sum(q3) - b*sum(q2)) / m)
        p(x)    = a*(x^b)
        error   = sum((Y - p.(X)) .^ 2)
        return p, error
    end
end