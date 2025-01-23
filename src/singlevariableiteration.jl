module SingleVariableIterations

export SingleVariableIteration
# export solve
export bisection
export fixed_point
export newton_raphson
export secant_method
export false_position

import ..LUSE_ENGR701_704_NumericalMethods: solve, newton_raphson, IntervalBoundsError

# KEY: [10.1145/3511528.3511535]
using Symbolics: @variables, Differential, simplify, expand_derivatives, build_function

# Ch. 2 (p. 47)
"""
    SingleVariableIteration(f, a, b, n, tol)

Given `f`(p) such that p ∈ [`a`, `b`], find the root of a single-variable equation in so many iterations, `n` within tolerance, `tol`.
"""
struct SingleVariableIteration{T<:Real}
    f   ::Function
    a   ::T
    b   ::T
    n   ::Int64
    tol ::T
end

"""
    maximum_slope(SVI::SingleVariableIteration)

Find the greatest value for first derivative of function.
"""
function maximum_slope(SVI::SingleVariableIteration)::AbstractFloat
    @variables x
    Dx  = Differential(x)
    df  = simplify(expand_derivatives(Dx(SVI.f(x))); expand=true)
    df  = build_function(df, x, expression=Val{false})
    return maximum(abs.(df.(range(SVI.a, SVI.b, length=1000))))
end

"""
    maximum_iterations(obj, method[, p0=0])

Find greatest integer for maximum iterations within tolerance.

# Notes
Acceptable values for `method` ∈ {`:bisection`, `:fixed_point`, `:newton_raphson`, `:secant_method`, `:false_position`}.
Initial guess, `p0` for function solution is not needed for `:bisection` method.
"""
function maximum_iterations(SVI::SingleVariableIteration, method::Symbol, p0::Float64 = 0.;
        k::Float64=NaN)::Int64
    if method == :bisection
        return ceil(-log(SVI.tol / (SVI.b - SVI.a))
            / log(2))
    elseif method ∈ (:fixed_point, :newton_raphson, :secant_method, :false_position)
        return ceil(-log(SVI.tol / max(p0 - SVI.a, SVI.b - p0))
            / log(k))
    else
        error("Method must be: `:bisection`, `:fixed_point`, `:newton_raphson`, `:secant_method`, or `:false_position`.")
    end
end

## 2.1 (p. 48)
"""
    solve(SVI::SingleVariableIteration[; method=:bisection, p0, p1, df])

Attempt to find where f(p) = 0 according to `method` ∈ {`:bisection` (default), `:fixed_point`, `:newton_raphson`, `:secant_method`, `:false_position`}.
Each `method` also has a convenience function.

# Notes
Convergence Rates:
    `:newton_raphson` > `:secant_method` > `:false_position` > `:fixed_point` > `:bisection`

`:bisection` is default because will always converge to a solution but is not the quickest.
"""
function solve(SVI::SingleVariableIteration;
    method  ::Symbol                    = :bisection,
    p0      ::T                         = 0.,
    p1      ::T                         = 0.,
    df      ::Union{Nothing, Function}  = nothing
)::Float64 where {T<:Float64}
    f, a, b = SVI.f, SVI.a, SVI.b
    if method ∈ (:bisection, :secant_method, :false_position)
        # check for opposite signs
        if (method == :bisection ? f(a)*f(b) : f(p0)*f(p1)) < 0
            k, N        = 1, SVI.n
            g, r        = Vector{Float64}(undef, N), Vector{Float64}(undef, N)
            g[k], r[k]  = f(a), 1.
            # exit by whichever condition is `true` first
            while r[k] >= SVI.tol && k < N
                if method == :bisection
                    x       = (b - a) / 2.
                    p       = a + x                         # new value, p
                    f(a)*f(p) > 0. ? (a = p) : (b = p)      # adjust next bounds
                    g[k+1], r[k+1] = p, abs(x)              # error of new value, p
                elseif method ∈ (:secant_method, :false_position)
                    q₀, q₁  = f(p0), f(p1)
                    p       = p1 - q₁*(p1 - p0)/(q₁ - q₀)   # new value, p
                    g[k+1], r[k+1] = p, abs(p - p0)         # error of new value
                    if method == :secant_method
                        p0      = p1
                    elseif method == :false_position
                        if f(p)*q₁ < 0.
                            p0 = p1                         # adjust next bounds
                        end
                    end
                    p1      = p
                end
                k += 1                                  # iterate to k + 1
            end
            return (k <= N && isassigned(g, k)) ? g[k] : NaN
        else # abort if f(a) is not opposite f(b)
            throw(IntervalBoundsError(SVI.f, SVI.a, SVI.b))
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
        k, N = 1, SVI.n
        g, r = Vector{Float64}(undef, N), Vector{Float64}(undef, N)
        g[k] = f(if method == :fixed_point
            (a + b) / 2.
        elseif method == :newton_raphson
            a
        end)
        r[k] = 1.
        # exit by whichever condition is `true` first
        while r[k] >= SVI.tol && k < N
            p = if method == :fixed_point
                f(p0)                           # new value, p
            elseif method == :newton_raphson
                p0 - (f(p0) / df(p0))
            end
            g[k+1], r[k+1] = p, if method == :fixed_point
                abs((p - p0) / p0)              # error of new value, p
            elseif method == :newton_raphson
                abs(p - p0)
            end
            p0 = p; k += 1                      # iterate to k + 1
        end
        return (k <= N && isassigned(g, k)) ? g[k] : NaN
    else
        throw(ArgumentError("The prescribed method must be one of the following selections: {`:bisection`, `:fixed_point`, `:newton_raphson`, `:secant_method`, `:false_position`}."))
    end
end

"""
    bisection(SVI::SingleVariableIteration)

Root-finding method: f(x) = 0.
Search for solution by halving the bounds such that `a` and `b` initially yield opposite signs in function.

# Notes
Relying on the **Intermediate Value Theorem** (IVT), this is a bracketed, root-finding method.
This method is rather slow to converge but will always converge to a solution; therefore, is a good starter method.
"""
bisection(SVI::SingleVariableIteration) = solve(SVI; method=:bisection)

## 2.2 (p. 55)
"""
    fixed_point(SVI::SingleVariableIteration, p0)

Attempt root-finding method with initial guess, `p0` in [a, b] by solving the equation g(p) = p via f(p) - p = 0.

**Use function with lowest slope!**

_Not root-bracketed._

# Notes
Theorem:
1) Existence of a fixed-point:
    If g ∈ C[a,b] and g(x) ∈ C[a, b] for all x ∈ [a, b], then function, g has a fixed point, p ∈ [a, b].
2) Uniqueness of a fixed point:
    If g'(x) exists on [a, b] and a positive constant, k < 1 exist with {|g'(x)| ≤ k | x ∈ (a, b)}, then there is exactly one fixed-point, p ∈ [a, b].
Converges by ``\\mathcal{O}(\\text{linear})`` if g'(p) ≠ 0, and ``\\mathcal{O}(\\text{quadratic})`` if g'(p) = 0 and g''(p) < M, where M = g''(ξ) that is the error function.
"""
fixed_point(SVI::SingleVariableIteration, p0::Float64) = solve(SVI; method=:fixed_point, p0=p0)

## 2.3 (p. 66)
"""
    newton_raphson(SVI::SingleVariableIteration, p0[; df=nothing])

Attempt root-finding method with initial guess, `p0` in [a, b] by solving the equation g(p) = p via f(p) - p = 0.
`df` will be the first derivative of function if not given.

**Use function with lowest slope!**

_`df`(x) ≠ 0_

# Notes
Quickest convergence rate, but not root-bracketed and has trouble with symmetric functions!
Initial guess, `p0` must be close to real solution; else, will converge to different root or oscillate (if symmetric).
This method can be viewed as fixed-point iteration.

Technique based on first Taylor polynomial expansion of f about ``p_{0}`` (that is `p0`) and evaluated at x = p.
|p - ``p_{0}``| is assumed small; therefore, 2ⁿᵈ-order Taylor term, the error, is small.

See `fixed_point()` for theorem.
"""
newton_raphson(SVI::SingleVariableIteration, p0::Float64;
    df::Union{Nothing, Function}=nothing
) = solve(SVI; method=:newton_raphson, p0=p0, df=df)

"""
    secant_method(SVI::SingleVariableIteration, p0, p1)

Attempt root-finding method with initial guesses such that `p0` and `p1` in [a, b] yield opposite signs in function.

**Use function with lowest slope!**

_Not root-bracketed._

# Notes
Method is less computationally expensive than `newton_raphson()` but may converge at slower rate by circumventing need to calculate derivative.

See `fixed_point()` for theorem.
"""
secant_method(SVI::SingleVariableIteration, p0::Float64, p1::Float64
) = solve(SVI; method=:secant_method, p0=p0, p1=p1)

"""
    false_position(SVI::SingleVariableIteration, p0, p1)

Attempt root-finding method with initial guesses such that `p0` and `p1` in [a, b] yield opposite signs in function.

**Use function with lowest slope!**

# Notes
Similar to, but slower to converge than, the `secant_method()` by including a test to ensure solution is root-bracketed.

See `fixed_point()` for theorem.
"""
false_position(SVI::SingleVariableIteration, p0::Float64, p1::Float64
) = solve(SVI; method=:false_position, p0=p0, p1=p1)

end
