module Integrations

export integrate

## 4.3 (p. 191)
"""
    integrate(f, x      [; rule=:trapezoidal, tol=10^-3])
    integrate(f, a, b, h[; rule=:trapezoidal, tol=10^-3])
    integrate(f, a, b, n[; rule=:trapezoidal, tol=10^-3])

Find the definite integral by some numerical quadrature.

# Notes
`f` may be a function or vector.
The domain may be defined with a vector, `x` or on the interval [`a`, `b`] either by number of sub-intervals, `n` or step-size, `h`.
`rule` accepts {`:trapezoidal` (default), `:midpoint`, `:simpson13`, `:simpson38`, `:simpsonN`}.
Dataset may contain unevenly spaces points.

# References
https://en.wikipedia.org/wiki/Simpson%27s_rule
"""
function integrate(f::Union{AbstractVector{T}, Function}, x::AbstractVector{T};
        rule::Symbol= :trapezoidal, tol::T= 10^-3)::Float64 where {T<:Real}
    is_function = isa(f, Function)
    @inbounds a, b, n     = x[begin], x[end], length(x) - 1
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
        @inbounds F           = integrate(is_function ? f : f[n:end], x[n:end], rule=:trapezoidal)
        @inbounds x           = x[begin:n]
        @inbounds a, b, n     = x[begin], x[end], length(x) - 1
    elseif rule == :simpson38 && n % 3 != 0
        m           = n - (n % 3 - 1)
        @inbounds F           = integrate(is_function ? f : f[m:end], x[m:end], rule=:trapezoidal)
        @inbounds x           = x[begin:m]
        @inbounds a, b, n     = x[begin], x[end], length(x) - 1
    elseif rule == :simpsonN && isodd(n)
        @inbounds hn2, hn1    = x[n] - x[n - 1], x[n + 1] - x[n]
        @inbounds fn2, fn1, fn = is_function ? f.(x[n - 1:n + 1]) : f[n - 1:n + 1]
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
            @inbounds z += (is_function ? f(x[j]) : f[j])
        end
        F          += if is_function
            h/2*(f(a) + 2z + f(b))
        else
            @inbounds h/2*(f[begin] + 2z + f[end])
        end
    elseif rule == :simpson13
        h           = (b - a) / n
        z1          = 0.
        for j ∈ 2:1:(n ÷ 2)
            @inbounds z1 += (is_function ? f(x[2j - 1]) : f[2j - 1])
        end
        z2          = 0.
        for j ∈ 1:1:(n ÷ 2)
            @inbounds z2 += is_function ? f(x[2j]) : f[2j]
        end
        F          += if is_function
            h/3*(f(a) + 2z1 + 4z2 + f(b))
        else
            @inbounds h/3*(f[begin] + 2z1 + 4z2 + f[end])
        end
    elseif rule == :simpson38
        h           = (b - a) / n
        z1          = 0.
        for j ∈ 2:1:n
            if j % 3 != 0
                @inbounds z1 += (is_function ? f(x[j]) : f[j])
            end
        end
        z3          = 0.
        for j ∈ 1:1:(n ÷ 3)
            @inbounds z3 += (is_function ? f(x[3j]) : f[3j])
        end
        F          += if is_function
            3h/8*(f(a) + 3z1 + 2z3 + f(b))
        else
            @inbounds 3h/8*(f[begin] + 3z1 + 2z3 + f[end])
        end
    elseif rule == :simpsonN
        h           = (b - a) / n
        for j ∈ 0:1:(n ÷ 2) - 1
            @inbounds h2j0, h2j1 = x[2j + 2] - x[2j + 1], x[2j + 3] - x[2j + 2]
            @inbounds f2j0, f2j1, f2j2 = is_function ? f.(x[2j + 1:2j + 3]) : f[2j + 1:2j + 3]
            F += 6 \ (h2j0 + h2j1) * (
                (2 - h2j1 / h2j0) * f2j0
                    + ((h2j0 * h2j1) \ (h2j0 + h2j1)^2) * f2j1
                    + (2 - h2j0 / h2j1) * f2j2)
        end
    elseif rule == :midpoint
        h           = (b - a) / (n + 2)
        z           = 0.
        for j ∈ 1:1:(n ÷ 2)
            @inbounds z += (is_function ? f(x[2j]) : f[2j])
        end
        F          += 2h*z
    else
        throw(ArgumentError("`rule` must be one of {`:trapezoidal`, `:simpson13`, `:simpson38`, `:simpsonN`, `:midpoint`}."))
    end
    return F
end

@inline function integrate(f::Union{AbstractVector{T}, Function}, a::T, b::T, h::T;
        rule::Symbol=:trapezoidal, tol::T=10^-3) where {T<:Real}
    return integrate(f, float.(a:h:b), rule=rule, tol=tol)
end

@inline function integrate(f::Union{AbstractVector{T}, Function}, a::T, b::T, n::Int64;
        rule::Symbol=:trapezoidal, tol::T=10^-3) where {T<:Real}
    return if rule == :midpoint
        integrate(f, float.(a:(b - a)/(n + 2):b), rule=rule, tol=tol)
    else
        integrate(f, float.(a:(b - a)/n:b), rule=rule, tol=tol)
    end
end

end
