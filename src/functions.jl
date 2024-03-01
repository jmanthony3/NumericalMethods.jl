using IntervalArithmetic, IntervalRootFinding
using Polynomials
using Printf
using StaticArrays
using Statistics
using Symbolics

# Ch. 2 (p. 47)
"""
    SingleVariableIteration(f, a, b, n, tol)

Given `f`(p) such that p ∈ [`a`, `b`], find the root of a single-variable equation in so many iterations, `n` within tolerance, `tol`.
"""
struct SingleVariableIteration
    f   ::Function
    a   ::Real
    b   ::Real
    n   ::Integer
    tol ::Real
end

"""
    maximum_slope(svi::SingleVariableIteration)

Find the greatest value for first derivative of function.
"""
function maximum_slope(svi::SingleVariableIteration)::AbstractFloat
    @variables x
    Dx  = Differential(x)
    df  = simplify(expand_derivatives(Dx(svi.f(x))); expand=true)
    df  = build_function(df, x, expression=Val{false})
    return maximum(abs.(df.(range(svi.a, svi.b, length=1000))))
end

"""
    maximum_iterations(obj, method[, p0=0])

Find greatest integer for maximum iterations within tolerance.

# Notes
Acceptable values for `method` ∈ {`:bisection`, `:fixed_point`, `:newton_raphson`, `:secant_method`, `:false_position`}.
Initial guess, `p0` for function solution is not needed for `:bisection` method.
"""
function maximum_iterations(
    svi     ::SingleVariableIteration,
    method  ::Symbol,
    p0      ::Real      = 0.;
    k       ::Real      = NaN
)::Integer
    if method == :bisection
        return ceil(-log(svi.tol / (svi.b - svi.a))
            / log(2))
    elseif method ∈ (:fixed_point, :newton_raphson, :secant_method, :false_position)
        return ceil(-log(svi.tol / max(p0 - svi.a, svi.b - p0))
            / log(k))
    else
        error("Method must be: `:bisection`, `:fixed_point`, `:newton_raphson`, `:secant_method`, or `:false_position`.")
    end
    # logging.info(f"With the inputs, I will terminate the technique after so many iterations, N = {max_iter}")
end

## 2.1 (p. 48)
"""
    solve(svi::SingleVariableIteration[; method=:bisection, p0, p1, df])

Attempt to find where f(p) = 0 according to `method` ∈ {`:bisection` (default), `:fixed_point`, `:newton_raphson`, `:secant_method`, `:false_position`}.
Each `method` also has a convenience function.

# Notes
Convergence Rates:
    `:newton_raphson` > `:secant_method` > `:false_position` > `:fixed_point` > `:bisection`

`:bisection` is default because will always converge to a solution but is not the quickest.
"""
function solve(svi::SingleVariableIteration;
    method  ::Symbol                    = :bisection,
    p0      ::T                         = 0.,
    p1      ::T                         = 0.,
    df      ::Union{Nothing, Function}  = nothing
)::Float64 where {T<:Float64}
    f, a, b = svi.f, svi.a, svi.b
    if method ∈ (:bisection, :secant_method, :false_position)
        # check for opposite signs
        if (method == :bisection ? f(a)*f(b) : f(p0)*f(p1)) < 0
            k, N        = 1, svi.n
            # g, r        = zeros(MVector{N}), zeros(MVector{N})
            g, r        = Vector{Float64}(undef, N), Vector{Float64}(undef, N)
            g[k], r[k]  = f(a), 1.
            # exit by whichever condition is `true` first
            while r[k] >= svi.tol && k < N
                if method == :bisection
                    x       = (b - a) / 2.
                    p       = a + x                         # new value, p
                    f(a)*f(p) > 0 ? (a = p) : (b = p)       # adjust next bounds
                    g[k+1], r[k+1] = p, abs(x)              # error of new value, p
                elseif method ∈ (:secant_method, :false_position)
                    q₀, q₁  = f(p0), f(p1)
                    p       = p1 - q₁*(p1 - p0)/(q₁ - q₀)   # new value, p
                    g[k+1], r[k+1] = p, abs(p - p0)         # error of new value
                    if method == :secant_method
                        p0      = p1
                    elseif method == :false_position
                        f(p)*q₁ < 0 ? p0 = p1 : nothing     # adjust next bounds
                    end
                    p1      = p
                end
                k += 1                                  # iterate to k + 1
            end
            return k <= N ? g[k] : NaN
        else # abort if f(a) is not opposite f(b)
            if method == :bisection
                error(@sprintf("Interval bounds must yield opposite signs in function, f := [f(a = %1.4f) = %1.4f, f(b = %1.4f) = %1.4f]",
                    a, b, f(a), f(b)))
            elseif method ∈ (:secant_method, :false_position)
                error(@sprintf("Interval bounds must yield opposite signs in function, f := [f(p₀ = %1.4f) = %1.4f, f(p₁ = %1.4f) = %1.4f]",
                    p0, p1, f(p0), f(p1)))
            end
        end
    elseif method ∈ (:fixed_point, :newton_raphson)
        # determine form of derivative
        if method == :newton_raphson && isnothing(df)
            @variables x
            Dx  = Differential(x)
            df  = simplify(expand_derivatives(Dx(f(x))); expand=true)
            df  = build_function(df, x, expression=Val{false})
        end
        # initialize
        k, N = 1, svi.n
        # g, r = zeros(MVector{N}), zeros(MVector{N})
        g, r = Vector{Float64}(undef, N), Vector{Float64}(undef, N)
        g[k] = f(if method == :fixed_point
            (a + b) / 2.
        elseif method == :newton_raphson
            a
        end)
        r[k] = 1.
        # exit by whichever condition is `true` first
        while r[k] >= svi.tol && k < N
            p = if method == :fixed_point
                f(p0)                           # new value, p
            elseif method == :newton_raphson
                p0 - (f(p0) / df(p0))
            end
            g[k+1], r[k+1] = p, if method == :fixed_point
                abs((p-p0)/p0)                  # error of new value, p
            elseif method == :newton_raphson
                abs(p - p0)
            end
            p0 = p; k += 1                      # iterate to k + 1
        end
        return k <= N ? g[k] : NaN
    end
end

"""
    bisection(svi::SingleVariableIteration)

Root-finding method: f(x) = 0.
Search for solution by halving the bounds such that `a` and `b` initially yield opposite signs in function.

# Notes
Relying on the Intermediate Value Theorem (IVT), this is a bracketed, root-finding method.
This method is rather slow to converge but will always converge to a solution; therefore, is a good starter method.
"""
bisection(svi::SingleVariableIteration) = solve(svi; method=:bisection)

## 2.2 (p. 55)
"""
    fixed_point(svi::SingleVariableIteration, p0)

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
fixed_point(svi::SingleVariableIteration, p0::Float64) = solve(svi; method=:fixed_point, p0=p0)

## 2.3 (p. 66)
"""
    newton_raphson(svi::SingleVariableIteration, p0[; df=nothing])

Attempt root-finding method with initial guess, `p0` in [a, b] by solving the equation g(p) = p via f(p) - p = 0. `df` will be the first derivative of function if not given.

**Use function with lowest slope!**

_`df`(x) ≠ 0_

# Notes
Quickest convergence rate, but not root-bracketed and has trouble with symmetric functions!
Initial guess, `p0` must be close to real solution; else, will converge to different root or oscillate (if symmetric).
This method can be viewed as fixed-point iteration.

Technique based on first Taylor polynomial expansion of f about p_{0} (that is `p0`) and evaluated at x = p. |p - p_{0}| is assumed small; therefore, 2^{\text{nd}}-order Taylor term, the error, is small.

See `fixed_point()` for theorem.
"""
newton_raphson(svi::SingleVariableIteration, p0::Float64;
    df::Union{Nothing, Function}=nothing
) = solve(svi; method=:newton_raphson, p0=p0, df=df)

"""
    secant_method(svi::SingleVariableIteration, p0, p1)

Attempt root-finding method with initial guesses such that `p0` and `p1` in [a, b] yield opposite signs in function.

**Use function with lowest slope!**

_Not root-bracketed._

# Notes
Method is less computationally expensive than `newton_raphson()` but may converge at slower rate by circumventing need to calculate derivative.

See `fixed_point()` for theorem.
"""
secant_method(svi::SingleVariableIteration, p0::Float64, p1::Float64
) = solve(svi; method=:secant_method, p0=p0, p1=p1)

"""
    false_position(svi::SingleVariableIteration, p0, p1)

Attempt root-finding method with initial guesses such that `p0` and `p1` in [a, b] yield opposite signs in function.

**Use function with lowest slope!**

# Notes
Similar to, but slower to converge than, the `secant_method()` by including a test to ensure solution is root-bracketed.

See `fixed_point()` for theorem.
"""
false_position(svi::SingleVariableIteration, p0::Float64, p1::Float64
) = solve(svi; method=:false_position, p0=p0, p1=p1)

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
    n::Union{Integer, Nothing}  = nothing
)::Tuple{Function, Vector{Float64}} where {T<:Vector{Float64}}
    @variables t
    Dt = Differential(t)
    function coefficient(x, xₖ, t)
        num, den = Num[], Float64[]
        for xₗ ∈ x
            if isa(xₗ, Num) || xₗ != xₖ
                push!(num, (t - xₗ))
                push!(den, (xₖ - xₗ))
            end
        end
        return prod(num) / prod(den)
    end
    function error(n, ξ, t)
        s       = Num[]
        # ξ_error = zeros(MVector{n})
        ξ_error = Vector{Float64}(undef, n)
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
                R = Float64[]
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
    n       = (isnothing(n) ? length(x) - 1 : n)
    # terms   = Num[]
    terms   = Vector{Num}(undef, n + 1)
    # errors  = zeros(MVector{n + 1})
    errors  = Vector{Float64}(undef, n + 1)
    for k ∈ 1:1:n + 1
        push!(terms, f[k] * coefficient(x, x[k], t))
        # terms[k]    = f[k] * coefficient(x, x[k], t)
        errors[k]   = error(k, sum(terms[begin:k]), t)
    end
    p       = build_function(sum(terms), t, expression=Val{false})
    return p, errors
end

function lagrange(
    x::T,
    f::T;
    n::Union{Integer, Nothing}  = nothing
) where {T<:AbstractVector}
    lagrange(float(collect(x)), float(collect(f)); n=n)
end;

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
    α   ::Float64;
    dir ::Symbol        = :auto
)::Function where {T<:Vector{Float64}}
    if dir == :auto
        dir = (α <= median(x) ? :forward : :backward)
    end
    fterm(g, i ,j) = (g[i, j] - g[i - 1, j]) / (g[i, 1] - g[i - (j - 1), 1])
    m, n    = length(x), length(x) + 1
    # fxn     = MMatrix{m, n}(zeros((m, n)))
    fxn     = Matrix{Float64}(undef, m, n)
    # coeff   = zeros(MVector{m - 1})
    coeff   = Vector{Float64}(undef, m - 1)
    fxn[:, 1], fxn[:, 2] = x, f
    for j ∈ 2:1:m, i ∈ j:1:m
        fₖ = fterm(fxn, i, j)
        fxn[i, j + 1] = fₖ
        if dir == :forward && i == j
            coeff[j - 1] = fₖ
        elseif dir == :backward && i == m
            coeff[j - 1] = fₖ
        end
    end
    @variables t
    k, g, terms = (dir == :forward ? 1 : m), 0., 1.
    for c ∈ coeff
        terms  *= (t - x[k])
        g      += c * terms
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
)::Tuple{Vector{Float64}, AbstractArray{Function}} where {T<:Vector{Float64}}
    function _algorithm(g)
        Y               = g
        m, n            = length(Y), length(Y) - 1
        # STEP 1:   build list, h_i
        # H               = zeros(MVector{n})
        H               = Vector{Float64}(undef, n)
        for i ∈ 1:1:n
            H[i] = X[i+1] - X[i]
        end
        # STEP 2:   build list, alpha_i
        # A, ALPHA        = Y, zeros(MVector{m})
        A, ALPHA        = Y, Vector{Float64}(undef, m)
        # ALPHA[1]        = 3*(A[2] - A[1])/H[1] - 3*AP[1]
        # ALPHA[m]        = 3*AP[m] - 3*(A[m] - A[n])/H[n]
        for i ∈ 2:1:n
            ALPHA[i] = 3/H[i]*(A[i+1] - A[i]) - 3/H[i-1]*(A[i] - A[i-1])
        end
        # Algorithm 6.7 to solve tridiagonal
        # STEP 3:   define l, mu, and z first points
        # L, MU           = zeros(MVector{m}), zeros(MVector{m})
        L, MU           = Vector{Float64}(undef, m), Vector{Float64}(undef, m)
        # Z, C            = zeros(MVector{m}), zeros(MVector{m})
        Z, C            = Vector{Float64}(undef, m), Vector{Float64}(undef, m)
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
        # B, D            = zeros(MVector{n}), zeros(MVector{n})
        B, D            = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
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
    n, splines      = length(X) - 1, Vector{Function}(undef, n) # []
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
)::Tuple{Vector{Float64}, AbstractArray{Function}} where {T<:Vector{Float64}}
    function _algorithm(g, gp)
        Y, YP           = g, gp
        m, n            = length(Y), length(Y) - 1
        # STEP 1:   build list, h_i
        # H               = zeros(MVector{n})
        H               = Vector{Float64}(undef, n)
        for i ∈ 1:1:n
            H[i] = X[i+1] - X[i]
        end
        # STEP 2:   define alpha list endpoints
        # A, AP, ALPHA    = Y, YP, zeros(MVector{m})
        A, AP, ALPHA    = Y, YP, Vector{Float64}(undef, m)
        ALPHA[1]        = 3(A[2] - A[1])/H[1] - 3AP[1]
        ALPHA[m]        = 3AP[m] - 3(A[m] - A[n])/H[n]
        # STEP 3:   build list, alpha_i
        for i ∈ 2:1:n
            ALPHA[i] = 3/H[i]*(A[i+1] - A[i]) - 3/H[i-1]*(A[i] - A[i-1])
        end
        # Algorithm 6.7 to solve tridiagonal
        # STEP 4:   define l, mu, and z first points
        # L, MU           = zeros(MVector{m}), zeros(MVector{m})
        L, MU           = Vector{Float64}(undef, m), Vector{Float64}(undef, m)
        # Z, C            = zeros(MVector{m}), zeros(MVector{m})
        Z, C            = Vector{Float64}(undef, m), Vector{Float64}(undef, m)
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
        # B, D            = zeros(MVector{n}), zeros(MVector{n})
        B, D            = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
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
    n, splines      = length(X) - 1, Vector{Function}(undef, n) # []
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
)::Float64 where {T<:Vector{Float64}}
    @variables t
    Dt = Differential(t)
    function coefficient(x, xₖ, t)
        num, den = Num[], Float64[]
        for xₗ ∈ x
            if isa(xₗ, Num) || xₗ != xₖ
                push!(num, (t - xₗ))
                push!(den, (xₖ - xₗ))
            end
        end
        return prod(num) / prod(den)
    end
    n       = (isnothing(n) ? length(x) - 1 : n)
    gp      = 0.
    for k ∈ 1:1:n + 1
        Lₖ          = coefficient(x, x[k], t)
        Lₖp         = simplify(expand_derivatives(Dt(Lₖ)); expand=true)
        Lₖp_eval    = build_function(Lₖp, t, expression=Val{false})
        gp         += f[k]*Lₖp_eval(x[j])
    end
    return gp
end

function n1derivative(
    x::T,
    f::T,
    j::Integer;
    n::Union{Integer, Nothing}  = nothing
) where {T<:AbstractVector}
    n1derivative(float(collect(x)), float(collect(f)), j; n=n)
end;

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
    InitialValueProblem(f, a, b, h, α, N)

Structure of the boundary conditions to Initial-Value Problem (IVP) differential equation.

# Notes
Make sure the independent variable (e.g. time) is the first argument of `f`!
"""
struct InitialValueProblem
    f::Function
    a::Real
    h::Real
    α::Real
    N::Integer
end

"""
    solve(ivp::InitialValueProblem[; method=:forward_euler, tol=10^-3])

Solve `ivp` according to `method` ∈ {`:forward_euler` (default), `:backward_euler`, `:improved_euler`, `:modified_euler`, `:runge_kutta`}.

# Notes
Each method has an equivalent convenience function.
E.g. `solve(ivp; method=:runge_kutta)` ≡ `runge_kutta(ivp)`.
"""
function solve(ivp::InitialValueProblem;
method::Symbol=:forward_euler, tol::Real=10^-3)::AbstractVector
    f       = ivp.f
    t, h, w = ivp.a, ivp.h, ivp.α
    ea, eb, λ = 1/2, 1/2, 1
    g       = zeros(MVector{ivp.N + 1})
    g[1]    = w
    for i ∈ 1:1:ivp.N
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
        t   = ivp.a + i*h
    end
    return g
end

forward_euler(ivp::InitialValueProblem;
tol=10^-3) = solve(ivp; tol=tol, method=:forward_euler)
backward_euler(ivp::InitialValueProblem;
tol=10^-3) = solve(ivp; tol=tol, method=:backward_euler)
improved_euler(ivp::InitialValueProblem;
tol=10^-3) = solve(ivp; tol=tol, method=:improved_euler)
modified_euler(ivp::InitialValueProblem;
tol=10^-3) = solve(ivp; tol=tol, method=:modified_euler)
runge_kutta(ivp::InitialValueProblem;
tol=10^-3) = solve(ivp; tol=tol, method=:runge_kutta)

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