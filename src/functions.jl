using IntervalArithmetic, IntervalRootFinding
using Polynomials
using StaticArrays
using Statistics
using Symbolics

# Ch. 3 (p. 103)
## 3.1 (p. 104)
"""
    lagrange()

Given a domain and range, construct a Lagrangian polynomial.
Polynomial will quickly begin to oscillate for larger datasets.
"""
function lagrange(
    x       ::T,
    f       ::T,
    degree  ::Union{Integer, Nothing}   = nothing
) where {T<:AbstractVector}
    @variables t
    Dt = Differential(t)
    function coefficient(xₖ, x)
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
            # println(("\tg", i, n, g))
            # println(("\tgp", i, n, gp))
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
            # println(("\tcoeffs", i, n, coeffs_sorted))
            # println(("\tpolynomial", i, n, Polynomial(coeffs_sorted)))
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
            # println(("\tξ", i, n, simplify(prod(ξ); expand=true)))
            for k ∈ 1:1:i
                ξ = simplify(expand_derivatives(Dt(ξ)); expand=true)
                # println(("\tξ", i, k, n, ξ))
            end
            # println(("\tξ / n!", i, n, ξ / (factorial(n))))
            Dξ = maximum(build_function(ξ, t; expression=Val{false}).(x[begin:i]) ./ (factorial(n)))
            # println((i, n, Dξ, abs(gx)))
            ξ_error[i] = Dξ * abs(gx)
            # xi_err = maximum(xi_error)
            # # println((i, n, xi_err))
        end
        return maximum(abs.(ξ_error))
    end
    degree  = (isnothing(degree) ? length(x) - 1 : degree)
    terms   = []
    errors  = zeros(MVector{degree + 1})
    for k ∈ 1:1:degree + 1
        # println(("degree", k, degree))
        # println(("term", k, term(domain[k], range[k], t)))
        push!(terms, f[k] * coefficient(x[k], t))
        # println(("terms", k, sum(terms)))
        errors[k]   = error(k, sum(terms), t)
    end
    # return nothing, nothing
    polynomial = build_function(sum(terms), t, expression=Val{false})
    return polynomial, errors
    # return polynomial
end

## 3.3 (p. 122)
"""
    newtondifference(x, f, α[; var::Symbol=:x, dir::Symbol=:auto])

Given a domain and range, construct some polynomial by Newton's Divided Difference.
`'forward'` or `'backward'` construction. Will be chosen automatically if not specified.

# Notes
Direction will be chosen if not specified.
Polynomials best made with even spacing in `domain`; although, this is not completely necessary.
"""
function newtondifference(
    x   ::T,
    f   ::T,
    α   ::Real;
    var ::Symbol        = :x,
    dir ::Symbol        = :auto
) where {T<:AbstractVector}
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
            # println((i, j))
            fₖ = fterm(fxn, i, j)
            fxn[i, j + 1] = fₖ
            if dir == :forward && i == j
                # println((i, j, fₖ))
                # push!(coeff, fₖ)
                coeff[j - 1] = fₖ
            elseif dir == :backward && i == m
                # println((i, j, fₖ))
                # push!(coeff, fₖ)
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
    natural(domain, f, function_derivative)

The bookend polynomials do not assume the slope entering and exiting the interval as the derivative at the respective endpoint.
"""
function natural(
    domain              ::T,
    f                   ::T
)::Tuple{AbstractVector, AbstractArray} where {T<:AbstractVector}
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
    g, X            = f, domain
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
    clamped(domain, f, function_derivative)

The bookend polynomials will have the same slope entering and exiting the interval as the derivative at the respective endpoint.
"""
function clamped(
    domain              ::T,
    f                   ::T,
    function_derivative ::T
)::Tuple{AbstractVector, AbstractArray} where {T<:AbstractVector}
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
    g, X, gp        = f, domain, function_derivative
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
function bezier(x, y, xguides, yguides)
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

function n1derivative(
    x       ::T,
    f       ::T,
    j       ::Integer;
    degree  ::Union{Integer, Nothing}   = nothing
)::AbstractFloat where {T<:AbstractVector}
    @variables t
    Dt = Differential(t)
    function coefficient(xₖ, x)
        num, den = [], []
        for xₗ ∈ x
            if isa(xₗ, Num) || xₗ != xₖ
                push!(num, (t - xₗ))
                push!(den, (xₖ - xₗ))
            end
        end
        return prod(num) / prod(den)
    end
    degree = (isnothing(degree) ? length(x) - 1 : degree)
    gp = 0.
    for k ∈ 1:1:degree + 1
        Lₖ          = coefficient(x[k], x)
        Lₖp         = simplify(expand_derivatives(Dt(Lₖ)); expand=true)
        Lₖp_eval    = build_function(Lₖp, t, expression=Val{false})
        gp         += f[k]*Lₖp_eval(x[j])
    end
    return gp
end

"""
    endpoint(x, y, h, point[, point_type="three"])

Find the derivative of a bookend point at either `:begin` or `:end` of dataset.
Acceptable values for `point_type` include {"three", "five"}.
"""
function endpoint(
    x           ::T,
    y           ::T,
    h           ::Real,
    point       ::Symbol;
    point_type  ::String    = "three"
)::AbstractFloat where {T<:AbstractVector}
    i = (point == :begin ? 1 : (point == :end ? length(x) : nothing))
    if point == :begin
        if point_type == "three"
            (-3y[i] + 4y[i+1] - y[i+2]) / 2h
        elseif point_type == "five"
            (-25y[i] + 48y[i+1]
                - 36y[i+2] + 16y[i+3]
                    - 3y[i+4]) / 12h
        else
            NaN
        end
    elseif point == :end
        if point_type == "three"
            -(-3y[i] + 4y[i-1] - y[i-2]) / 2h
        elseif point_type == "five"
            -(-25y[i] + 48y[i-1]
                - 36y[i-2] + 16y[i-3]
                    - 3y[i-4]) / 12h
        else
            NaN
        end
    else
        NaN
    end
end

"""
    midpoint(x, y, h, point[, point_type="three"])

Find the derivative of some point within a dataset.
Acceptable values for `point_type` include {"three", "five", "2nd_derivative"}.
"""
function midpoint(
    x           ::T,
    y           ::T,
    h           ::Real,
    point       ::Integer;
    point_type  ::String    = "three"
)::AbstractFloat where {T<:AbstractVector}
    i = point
    return if point_type == "three"
        (y[i + 1] - y[i - 1]) / 2h
    elseif point_type == "five"
        (y[i - 2] - 8y[i - 1]
            + 8y[i + 1] - y[i + 2]) / 12h
    elseif point_type == "2nd_derivative"
        (y[i - 1] - 2y[i] + y[i + 1]) / (h ^ 2.)
    else
        NaN
    end
end

## 4.3 (p. 191)
"""
    integrate(f[; rule=:trapezoidal, tol=10^-3])

# Notes
Find the definite integral by some composite numeric quadrature.
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

# Ch. 5 (p. 259)
"""
    ODE(f, a, b, h, α, N)

Structure of the boundary conditions to differential equation where `f` is the time derivative of the function to approximate.

# Notes
Make sure the independent variable is the first argument of `f`!
"""
struct ODE
    f   ::Function
    a   ::Real
    b   ::Real
    h   ::Real
    α   ::Real
    N   ::Integer
    # vars::AbstractVector{Num}
end

"""
    ivp(obj::ODE[, tol=10^-3; method=:forward_euler])

Solve `obj` according to `method` ∈ {`:forward_euler` (default), `:backward_euler`, `:improved_euler`, `:modified_euler`, `:runge_kutta`}.

# Notes
Each method has an equivalent convenience function.
E.g. `ivp(obj; method=:runge_kutta)` ≡ `runge_kutta(obj)`.
"""
function ivp(obj::ODE, tol=10^-3; method=:forward_euler)
    f       = obj.f
    a, b, h = obj.a, obj.b, obj.h
    # vars    = obj.vars
    t, w0   = a, obj.α
    ea, eb, λ = 1/2, 1/2, 1
    g       = zeros(MVector{obj.N + 1})
    g[1]    = w0
    # g       = [w0]
    for i ∈ 1:1:obj.N
        w = w0 + if method == :foward_euler
            h * f(t, w0)
        elseif method == :backward_euler
            h * f(t + h, w0 + h*f(t, w0))
        elseif method ∈ [:improved_euler, :modified_euler]
            h * (ea*f(t, w0) + eb*f(t + λ*h, w0 + λ*h*f(t, w0)))
        elseif method == :runge_kutta
            k1 = h * f(t,       w0)
            k2 = h * f(t + h/2, w0 + k1/2)
            k3 = h * f(t + h/2, w0 + k2/2)
            k4 = h * f(t + h,   w0 + k3)
            (k1 + 2k2 + 2k3 + k4) / 6
        end
        # push!(g, w)
        g[i + 1] = w
        # push!(increments, w - w0)
        w0  = w
        t   = a + i*h
    end
    return g
end

forward_euler(obj, tol=10^-3) = ivp(obj, tol; method=:forward_euler)
backward_euler(obj, tol=10^-3) = ivp(obj, tol; method=:backward_euler)
improved_euler(obj, tol=10^-3) = ivp(obj, tol; method=:improved_euler)
modified_euler(obj, tol=10^-3) = ivp(obj, tol; method=:modified_euler)
runge_kutta(obj, tol=10^-3) = ivp(obj, tol; method=:runge_kutta)

# Ch. 8 (p. 505)
## 8.1 (p. 506)
"""
    linearleastsquares(domain, f, degree::Integer)

Construct a polynomial of some degree while minimizing the least squares error.

# Notes
Least squares error := ``E = \\sum_{i=1}^{m}[y_{i} - P_{n}(x_{i})]^{2}``

Constructed polynomial of the form: ``P(x) = a_{n}x^{n} + a_{n - 1}x^{n - 1} + \\dots + a_{1}x + a_{0}``
"""
function linearleastsquares(
    domain  ::T,
    f       ::T,
    degree  ::Integer
) where {T<:AbstractVector}
    X, Y            = domain, f
    m               = length(X)
    A, b            = zeros((degree+1, degree+1)), zeros(degree+1)
    for i ∈ 0:1:degree
        for j ∈ 0:1:degree
            for k ∈ 1:1:m
                A[i+1,j+1] += X[k]^(i + j)
            end
        end
        for j ∈ 1:1:m
            b[i+1] += Y[j]*X[j]^i
        end
    end
    @variables t
    x, polynomial   = A\b, 0
    for i ∈ 0:1:degree
        polynomial += x[i+1]*t^i
    end
    polynomial      = build_function(polynomial, t, expression=Val{false})
    error           = sum((Y - polynomial.(X)).^2)
    return polynomial, error
end

"""
    linearleastsquares(domain, f, type::Symbol)

Given a domain and range, yield the coefficients for an equation and the equation of the form ``y = ax^{b}``.
"""
function linearleastsquares(
    domain  ::T,
    f       ::T,
    type    ::Symbol
) where {T<:AbstractVector}
    if type == :power
        X, Y            = domain, f
        m               = length(X)
        q1, q2, q3, q4  = [], [], [], []
        for i ∈ 1:1:m
            push!(q1, log(X[i])*log(Y[i]))
            push!(q2, log(X[i]))
            push!(q3, log(Y[i]))
            push!(q4, log(X[i])^2)
        end
        num             = m*sum(q1) - sum(q2)*sum(q3)
        den             = m*sum(q4) - (sum(q2))^2
        b               = num / den
        a               = exp((sum(q3) - b*sum(q2)) / m)
        expression(x)   = a*(x^b)
        error           = sum((Y - expression.(X)) .^ 2)
        return expression, error
    end
end