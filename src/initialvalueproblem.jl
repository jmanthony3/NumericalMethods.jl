module InitialValueProblems

export InitialValueProblem
# export solve
export forward_euler
export backward_euler
export improved_euler
export modified_euler
export runge_kutta

import ..NumericalMethods: solve

using StaticArrays

# Ch. 5 (p. 259)
"""
    InitialValueProblem(f, a, b, h, α, N)

Structure of the boundary conditions to **Initial-Value Problem** (IVP) differential equation.

# Notes
Make sure the independent variable (e.g. time) is the first argument of `f`!
"""
struct InitialValueProblem{T<:Real}
    f::Function
    a::T
    h::T
    α::T
    N::Int64
end

"""
    solve(ivp::InitialValueProblem[; method=:forward_euler, tol=10^-3])

Solve `ivp` according to `method` ∈ {`:forward_euler` (default), `:backward_euler`, `:improved_euler`, `:modified_euler`, `:runge_kutta`}.

# Notes
Each `method` has an equivalent convenience function.
E.g. `solve(ivp; method=:runge_kutta)` ≡ `runge_kutta(ivp)`.
"""
function solve(ivp::InitialValueProblem;
        method::Symbol=:forward_euler, tol::Real=10^-3)::Vector{Float64}
    f       = ivp.f
    t, h, w = ivp.a, ivp.h, ivp.α
    ea, eb, λ = 1/2, 1/2, 1
    g       = zeros(MVector{ivp.N + 1})
    @inbounds g[1]    = w
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
        @inbounds g[i + 1] = w
        t   = ivp.a + i*h
    end
    return g
end

forward_euler(ivp::InitialValueProblem;     tol=10^-3) = solve(ivp; tol=tol, method=:forward_euler)
backward_euler(ivp::InitialValueProblem;    tol=10^-3) = solve(ivp; tol=tol, method=:backward_euler)
improved_euler(ivp::InitialValueProblem;    tol=10^-3) = solve(ivp; tol=tol, method=:improved_euler)
modified_euler(ivp::InitialValueProblem;    tol=10^-3) = solve(ivp; tol=tol, method=:modified_euler)
runge_kutta(ivp::InitialValueProblem;       tol=10^-3) = solve(ivp; tol=tol, method=:runge_kutta)

end
