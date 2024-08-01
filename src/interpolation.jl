module Interpolations

export linearinterpolation
export lagrange
export newtondifference
export natural
export clamped
export bezier
export linearleastsquares

using IntervalArithmetic, IntervalRootFinding # KEY: [IntervalArithmetic.jl]
using Polynomials: Polynomial, roots
using StaticArrays
using Statistics: median
# KEY: [10.1145/3511528.3511535]
using Symbolics: Num, @variables, Differential, simplify, expand_derivatives, build_function, value, degree

"""
    linearinterpolation(x0, y0, x1, y1, x)

``y = y₀ + (x - x₀)*(y₁ - y₀)/(x₁ - x₀)``
"""
@inline function linearinterpolation(x0::T, y0::T, x1::T, y1::T, x::T)::Float64 where {T<:Float64}
    return y0 + (x - x0)*(y1 - y0)/(x1 - x0)
end

"""
    linearinterpolation(x, y, p)

Calls `linearinterpolation` with the first index in `x` less than and greater than `p`.
If for any `p` ∈ `x`, then the first occurrence in `y` is returned.
"""
function linearinterpolation(x::Vector{T}, y::Vector{T}, p::T)::Float64 where {T<:Float64}
    return if isempty(findall(x->x==p, x))
        i, j    = findlast(x .< p), findfirst(x .> p)
        x0, y0  = x[i], y[i]
        x1, y1  = x[j], y[j]
        linearinterpolation(x0, y0, x1, y1, p)
    else
        y[findfirst(x .== p)]
    end
end

# Ch. 3 (p. 103)
## 3.1 (p. 104)
function lagrange_coefficient(x::Vector{T}, xₖ::T, t::Num) where {T<:Float64}
    num, den = Num[], Float64[]
    for xₗ ∈ x
        if isa(xₗ, Num) || xₗ != xₖ
            push!(num, (t - xₗ))
            push!(den, (xₖ - xₗ))
        end
    end
    return prod(num) / prod(den)
end

"""
    lagrange(x, f[; n=nothing])

Given a domain, `x` and range, `f`, construct the `n`th Lagrangian polynomial.

# Notes
If `n=nothing`, then method will utilize entire dataset.
Polynomials will quickly oscillate for larger datasets.
"""
function lagrange(x::T, f::T;
        n::Union{Int64, Nothing}=nothing)::Tuple{Function, Vector{Float64}} where {T<:Vector{Float64}}
    @variables t
    Dt = Differential(t)
    function error(n, ξ, t)
        s       = Num[]
        # ξ_error = zeros(MVector{n})
        ξ_error = Vector{Float64}(undef, n)
        for i ∈ 1:1:n
            @inbounds push!(s, t - x[i])
            g               = simplify(prod(s); expand=true)
            g_eval          = build_function(g, t, expression=Val{false})
            gp              = simplify(expand_derivatives(Dt(g)); expand=true)
            gp_eval         = build_function(gp, t, expression=Val{false})
            coeffs_sorted   = if hasproperty(value(gp), :dict)
                coeffs_map      = value(gp).dict
                exps            = degree.(keys(coeffs_map))
                collect(values(coeffs_map))[sortperm(exps)]
            else
                []
            end
            if gp_eval(0.) != 0.
                pushfirst!(coeffs_sorted, gp_eval(0.))
            end
            g_prime_poly    = Polynomial(float.(coeffs_sorted))
            gx              = if i == 1
                g_eval(first(coeffs_sorted))
            elseif i == 2
                g_eval(first(roots(g_prime_poly)))
            else
                R = Float64[]
                for r ∈ roots(g_prime_poly)
                    if isreal(r)
                        push!(R, g_eval(r))
                    end
                end
                maximum(abs.(R))
            end
            for _ ∈ 1:1:i
                ξ = simplify(expand_derivatives(Dt(ξ)); expand=true)
            end
            dξ              = maximum(
                build_function(ξ, t; expression=Val{false}).(x[begin:i])
                    ./ (factorial(n)))
            @inbounds ξ_error[i]      = dξ * abs(gx)
        end
        return maximum(abs.(ξ_error))
    end
    n       = (isnothing(n) ? length(x) - 1 : n)
    terms   = Vector{Num}(undef, n + 1)
    errors  = Vector{Float64}(undef, n + 1)
    for k ∈ 1:1:n + 1
        @inbounds terms[k]    = f[k] * lagrange_coefficient(x, x[k], t)
        @inbounds errors[k]   = error(k, sum(terms[begin:k]), t)
    end
    p       = build_function(sum(terms), t, expression=Val{false})
    return p, errors
end

@inline function lagrange(x::Vector{T}, f::Vector{T};
        n::Union{Int64, Nothing}=nothing) where {T<:Float64}
    return lagrange(float(collect(x)), float(collect(f)); n=n)
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
function newtondifference(x::Vector{T}, f::Vector{T}, α::Float64;
        dir::Symbol=:auto)::Function where {T<:Float64}
    if dir == :auto
        dir = (α <= median(x) ? :forward : :backward)
    end
    fterm(g, i ,j) = (g[i, j] - g[i - 1, j]) / (g[i, 1] - g[i - (j - 1), 1])
    m, n    = length(x), length(x) + 1
    fxn     = Matrix{Float64}(undef, m, n)
    coeff   = Vector{Float64}(undef, m - 1)
    fxn[:, 1], fxn[:, 2] = x, f
    for j ∈ 2:1:m, i ∈ j:1:m
        fₖ = fterm(fxn, i, j)
        @inbounds fxn[i, j + 1] = fₖ
        if dir == :forward && i == j
            @inbounds coeff[j - 1] = fₖ
        elseif dir == :backward && i == m
            @inbounds coeff[j - 1] = fₖ
        end
    end
    @variables t
    k, g, terms = (dir == :forward ? 1 : m), 0., 1.
    for c ∈ coeff
        @inbounds terms  *= (t - x[k])
        g      += c * terms
        k      += (dir == :forward ? 1 : -1)
    end
    p = g + (dir == :forward ? first(f) : last(f))
    return build_function(p, t, expression=Val{false})
end

## 3.5 (p. 142)
"""
    natural(x, f)

The bookend polynomials do not assume the slope entering and exiting the interval as the derivative at the respective endpoint.
"""
function natural(x::Vector{T}, f::Vector{T}
        )::Tuple{Vector{Float64}, Vector{Function}} where {T<:Float64}
    function _algorithm(g)
        Y               = g
        m, n            = length(Y), length(Y) - 1
        # STEP 1:   build list, h_i
        H               = Vector{Float64}(undef, n)
        for i ∈ 1:1:n
            @inbounds H[i] = X[i+1] - X[i]
        end
        # STEP 2:   build list, alpha_i
        A, ALPHA        = Y, Vector{Float64}(undef, m)
        for i ∈ 2:1:n
            @inbounds ALPHA[i] = 3/H[i]*(A[i+1] - A[i]) - 3/H[i-1]*(A[i] - A[i-1])
        end
        # Algorithm 6.7 to solve tridiagonal
        # STEP 3:   define l, mu, and z first points
        L, MU           = Vector{Float64}(undef, m), Vector{Float64}(undef, m)
        Z, C            = Vector{Float64}(undef, m), Vector{Float64}(undef, m)
        @inbounds L[1], MU[1], Z[1]= 1., 0, 0.
        # STEP 4:   build lists l, mu, and z
        for i ∈ 2:1:n
            @inbounds L[i]  = 2(X[i+1] - X[i-1]) - H[i-1]*MU[i-1]
            @inbounds MU[i] = H[i] / L[i]
            @inbounds Z[i]  = (ALPHA[i] - H[i-1]*Z[i-1]) / L[i]
        end
        # STEP 5:   define l, z, and c endpoints
        @inbounds L[m], Z[m], C[m]= 1., 0., 0.
        # STEP 6:   build lists c, b, and d
        B, D            = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
        for i ∈ 0:1:n-1
            j    = n-i
            @inbounds C[j] = Z[j] - MU[j]*C[j+1]
            @inbounds B[j] = (A[j+1] - A[j])/H[j] - H[j]*(C[j+1] + 2C[j])/3
            @inbounds D[j] = (C[j+1] - C[j]) / 3H[j]
        end
        return Y, A, B, C, D
    end
    g, X            = f, x
    Y, A, B, C, D   = _algorithm(g)
    n               = length(X) - 1
    splines         = Vector{Function}(undef, n) # []
    for j ∈ 1:1:n
        @inbounds xj, aj, bj, cj, dj = X[j], A[j], B[j], C[j], D[j]
        sj(x) = aj + bj*(x - xj) + cj*(x - xj)^2 + dj*(x - xj)^3
        @inbounds splines[j] = sj
    end
    return Y, splines
end

"""
    clamped(x, f, fp)

The bookend polynomials will have the same slope entering and exiting the interval as the derivative at the respective endpoint.
"""
function clamped(x::Vector{T}, f::Vector{T}, fp::Vector{T}
        )::Tuple{Vector{Float64}, Vector{Function}} where {T<:Float64}
    function _algorithm(g, gp)
        Y, YP           = g, gp
        m, n            = length(Y), length(Y) - 1
        # STEP 1:   build list, h_i
        H               = Vector{Float64}(undef, n)
        for i ∈ 1:1:n
            @inbounds H[i] = X[i+1] - X[i]
        end
        # STEP 2:   define alpha list endpoints
        A, AP, ALPHA    = Y, YP, Vector{Float64}(undef, m)
        @inbounds ALPHA[1]        = 3(A[2] - first(A))/first(H) - first(3AP)
        @inbounds ALPHA[m]        = 3AP[m] - 3(A[m] - A[n])/H[n]
        # STEP 3:   build list, alpha_i
        for i ∈ 2:1:n
            @inbounds ALPHA[i] = 3/H[i]*(A[i+1] - A[i]) - 3/H[i-1]*(A[i] - A[i-1])
        end
        # Algorithm 6.7 to solve tridiagonal
        # STEP 4:   define l, mu, and z first points
        L, MU           = Vector{Float64}(undef, m), Vector{Float64}(undef, m)
        Z, C            = Vector{Float64}(undef, m), Vector{Float64}(undef, m)
        @inbounds L[1], MU[1]     = 2first(H), 0.5
        @inbounds Z[1]            = first(ALPHA) / first(L)
        # STEP 5:   build lists l, mu, and z
        for i ∈ 2:1:n
            @inbounds L[i]  = 2(X[i+1] - X[i-1]) - H[i-1]*MU[i-1]
            @inbounds MU[i] = H[i]/L[i]
            @inbounds Z[i]  = (ALPHA[i] - H[i-1]*Z[i-1])/L[i]
        end
        # STEP 6:   define l, z, and c endpoints
        @inbounds L[m]            = H[n] * (2 - MU[n])
        @inbounds Z[m]            = (ALPHA[m] - H[n]*Z[n]) / L[m]
        @inbounds C[m]            = Z[m]
        # STEP 7:   build lists c, b, and d
        B, D            = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
        for i ∈ 0:1:n-1
            j    = n-i
            @inbounds C[j] = Z[j] - MU[j]*C[j+1]
            @inbounds B[j] = (A[j+1] - A[j])/H[j] - H[j]*(C[j+1] + 2*C[j])/3
            @inbounds D[j] = (C[j+1] - C[j]) / 3H[j]
        end
        return Y, A, B, C, D
    end
    g, X, gp        = f, x, fp
    Y, A, B, C, D   = _algorithm(g, gp)
    n               = length(X) - 1
    splines         = Vector{Function}(undef, n)
    for j ∈ 1:1:n
        @inbounds xj, aj, bj, cj, dj = X[j], A[j], B[j], C[j], D[j]
        sj(x) = aj + bj*(x - xj) + cj*(x - xj)^2 + dj*(x - xj)^3
        @inbounds splines[j] = sj
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
function bezier(x::T, y::T, xguides::T, yguides::T)::Vector{Function} where {T<:Vector{Float64}}
    n, curves = length(x) - 1, []
    for i ∈ 1:1:n
        @inbounds a = (x[i],
            3(xguides[i] - x[i]),
            3(x[i] + xguides[i + 1] - 2xguides[i]),
            x[i + 1] - x[i] + 3xguides[i] - 3xguides[i + 1])
        @inbounds b = (y[i],
            3(yguides[i] - y[i]),
            3(y[i] + yguides[i + 1] - 2yguides[i]),
            y[i + 1] - y[i] + 3yguides[i] - 3yguides[i + 1])
        xi(t) = sum(a .* (1., t, t^2, t^3))
        yi(t) = sum(b .* (1., t, t^2, t^3))
        push!(curves, (x = xi, y = yi))
    end
    return curves
end

# Ch. 8 (p. 505)
## 8.1 (p. 506)
"""
    linearleastsquares(x, f, n::Int64)

Construct a polynomial of degree, `n` while minimizing the least squares error.

# Notes
Least squares error := ``E = \\sum_{i=1}^{m}[y_{i} - P_{n}(x_{i})]^{2}``

Constructed polynomial of the form: ``P(x) = a_{n}x^{n} + a_{n - 1}x^{n - 1} + \\dots + a_{1}x + a_{0}``
"""
function linearleastsquares(x::T, f::T, n::Int64)::Tuple{Function, Float64} where {T<:Vector{Float64}}
    X, Y    = x, f
    m       = length(X)
    A       = MMatrix{n+1, n+1}(zeros((n+1, n+1)))
    b       = zeros(MVector{n+1})
    for i ∈ 0:1:n
        for j ∈ 0:1:n
            for k ∈ 1:1:m
                @inbounds A[i+1,j+1] += X[k]^(i + j)
            end
        end
        for j ∈ 1:1:m
            @inbounds b[i+1] += Y[j]*X[j]^i
        end
    end
    @variables t
    x, p    = A\b, 0.
    for i ∈ 0:1:n
        @inbounds p  += x[i+1]*t^i
    end
    p       = build_function(p, t, expression=Val{false})
    error   = sum((Y - p.(X)) .^ 2)
    return p, error
end

"""
    linearleastsquares(x, f, type::Symbol)

Given a domain and range, yield the coefficients for an equation and the equation of the form ``y = ax^{b}``.
"""
function linearleastsquares(x::T, f::T, type::Symbol)::Tuple{Function, Float64} where {T<:Vector{Float64}}
    if type == :power
        X, Y    = x, f
        m       = length(X)
        q1      = zeros(MVector{m})
        q2      = zeros(MVector{m})
        q3      = zeros(MVector{m})
        q4      = zeros(MVector{m})
        for i ∈ 1:1:m
            @inbounds q1[i] = log(X[i])*log(Y[i])
            @inbounds q2[i] = log(X[i])
            @inbounds q3[i] = log(Y[i])
            @inbounds q4[i] = log(X[i])^2
        end
        num     = m*sum(q1) - sum(q2)*sum(q3)
        den     = m*sum(q4) - sum(q2)^2
        b       = num / den
        a       = exp((sum(q3) - b*sum(q2)) / m)
        p(x)    = a*(x^b)
        error   = sum((Y - p.(X)) .^ 2)
        return p, (isnothing(error) ? 0. : error)
    end
end

end
