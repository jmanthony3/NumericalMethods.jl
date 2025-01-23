module BoundaryValueProblems

export BoundaryValueProblem
# export solve
export linear_shooting_method
export finite_difference_method

import ..LUSE_ENGR701_704_NumericalMethods: solve, SystemOfEquation

using LinearAlgebra

# Ch. 11 (p. 685)
"""
    BoundaryValueProblem(f, a, b, h, α, β, N)

Structure of the boundary conditions to **Boundary-Value Problem** (BVP) differential equations.

# Notes
Make sure the independent variable (e.g. time) is the first argument of `f`!
"""
struct BoundaryValueProblem{T<:Real}
    f::Vector{Function}
    a::T
    b::T
    h::T
    α::T
    β::T
    N::Int64
end

## 11.1 (p. 686)
"""
    linear_shooting_method(BVP)

Solve a BVP differential equation with 2 IVP differential equations with in RK4 scheme.
"""
function linear_shooting_method(BVP)::Vector{NTuple{2, Float64}}
    u1, u2, v1, v2 = [BVP.α], [0.], [0.], [1.]
    p, q, r = BVP.f
    for (i, t) in enumerate(BVP.α:BVP.h:BVP.β-BVP.h)
        k11 = BVP.h*u2[i]
        k12 = BVP.h*(p(t)*u2[i] + q(t)*u1[i] + r(t))
        k21 = BVP.h*(u2[i] + k12/2)
        k22 = BVP.h*(p(t + BVP.h/2)*(u2[i] + k12/2) + q(t + BVP.h/2)*(u1[i] + k11/2) + r(t + BVP.h/2))
        k31 = BVP.h*(u2[i] + k22/2)
        k32 = BVP.h*(p(t + BVP.h/2)*(u2[i] + k22/2) + q(t + BVP.h/2)*(u1[i] + k21/2) + r(t + BVP.h/2))
        k41 = BVP.h*(u2[i] + k32)
        k42 = BVP.h*(p(t + BVP.h)*(u2[i] + k32) + q(t + BVP.h)*(u1[i] + k31) + r(t + BVP.h))
        push!(u1, u1[i] + (k11 + 2k21 + 2k31 + k41) / 6)
        push!(u2, u2[i] + (k12 + 2k22 + 2k32 + k42) / 6)
        ###############################
        k11 = BVP.h*v2[i]
        k12 = BVP.h*(p(t)*v2[i] + q(t)*v1[i])
        k21 = BVP.h*(v2[i] + k12/2)
        k22 = BVP.h*(p(t + BVP.h/2)*(v2[i] + k12/2) + q(t + BVP.h/2)*(v1[i] + k11/2))
        k31 = BVP.h*(v2[i] + k22/2)
        k32 = BVP.h*(p(t + BVP.h/2)*(v2[i] + k22/2) + q(t + BVP.h/2)*(v1[i] + k21/2))
        k41 = BVP.h*(v2[i] + k32)
        k42 = BVP.h*(p(t + BVP.h)*(v2[i] + k32) + q(t + BVP.h)*(v1[i] + k31))
        push!(v1, v1[i] + (k11 + 2k21 + 2k31 + k41) / 6)
        push!(v2, v2[i] + (k12 + 2k22 + 2k32 + k42) / 6)
    end
    w1, w2 = [BVP.α], [(BVP.β - last(u1)) / last(v1)]
    for i in 2:1:BVP.N+1
        push!(w1, u1[i] + first(w2)*v1[i])
        push!(w2, u2[i] + first(w2)*v2[i])
    end
    return [(y, yp) for (y, yp) in zip(w1, w2)]
end

## 11.3 (p. 700)
"""
    finite_difference_method(BVP[; method=:gauss_seidel, M=100, tol=10^-3])

Solve BVP differential equations with Dirichlet boundary conditions by 2 IVP differential equations according to `method` ∈ {`:jacobi`, `:gauss_seidel` (default), `:successive_relaxation`, `:newton_raphson`} within numerical iterations, `M` and tolerance, `tol`.

Uses a Taylor polynomial with a first-order and a second-order IVP equations.
*Converges ``\\mathcal{O}(h^{2})``*.
"""
function finite_difference_method(BVP::BoundaryValueProblem; method::Symbol=:gauss_seidel, M::Int64=100, tol::Real=10^-6)::Vector{Float64}
    p, q, r = BVP.f
    t = BVP.a + BVP.h
    ai = [2 + (BVP.h^2.) * q(t)]
    bi = [-1 + (BVP.h/2) * p(t)]
    ci = Float64[]
    di = [-(BVP.h^2.) * r(t) + (1 + (BVP.h / 2) * p(t)) * BVP.α]
    for t in BVP.a+2BVP.h:BVP.h:BVP.b-2BVP.h
        push!(ai, 2 + (BVP.h^2) * q(t))
        push!(bi, -1 + (BVP.h/2) * p(t))
        push!(ci, -1 - (BVP.h/2) * p(t))
        push!(di, -(BVP.h^2) * r(t))
    end
    t = BVP.b - BVP.h
    push!(ai, 2 + (BVP.h^2) * q(t))
    push!(ci, -1 - (BVP.h/2) * p(t))
    push!(di, -(BVP.h^2) * r(t) + (1 - (BVP.h / 2) * p(t)) * BVP.β)
    # MVI = MultiVariableIteration(diagm(-1 => ci, 0 => ai, 1 => bi), zeros(length(ai)), di, M, tol)
    # return solve(MVI; method=method)
    return [BVP.α, solve(SystemOfEquation(diagm(-1 => ci, 0 => ai, 1 => bi), di, M, tol))..., BVP.β]
end

# p(x) = -x\2
# q(x) = 2/(x^2)
# r(x) = sin(log(x))/x^2
# BVP = BoundaryValueProblem([p, q, r], 1., 2., 0.1, 1., 2., 10)
# for (y, yp) in (linear_shooting_method(BVP))
#     println((y, yp))
# end
# for x in finite_difference_method(BVP; tol=10^-6)
#     println(x)
# end

end
