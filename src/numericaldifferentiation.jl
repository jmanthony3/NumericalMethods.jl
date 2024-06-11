module Derivatives

export n1derivative
export endpoint
export midpoint

import ..Interpolations: lagrange_coefficient

using Symbolics

# Ch. 4 (p. 171)
## 4.1 (p. 172)
"""
    n1derivative(x, f, j[; n=nothing])

The general (n + 1)-point formula to approximate f' at point `j`.

# Notes
If `n = nothing`, then entire dataset used to construct `n`'th Lagrange coefficient.
"""
function n1derivative(x::Vector{T}, f::Vector{T}, j::Int64;
        n::Union{Int64, Nothing}=nothing)::Float64 where {T<:Float64}
    @variables t
    Dt = Differential(t)
    n       = (isnothing(n) ? length(x) - 1 : n)
    gp      = 0.
    for k ∈ 1:1:n + 1
        @inbounds Lₖ    = lagrange_coefficient(x, x[k], t)
        Lₖp             = simplify(expand_derivatives(Dt(Lₖ)); expand=true)
        Lₖp_eval        = build_function(Lₖp, t, expression=Val{false})
        @inbounds gp   += f[k]*Lₖp_eval(x[j])
    end
    return gp
end

@inline function n1derivative(x::AbstractVector{T}, f::AbstractVector{T}, j::Int64;
        n::Union{Int64, Nothing}=nothing) where {T<:Float64}
    return n1derivative(float(collect(x)), float(collect(f)), j; n=n)
end

"""
    endpoint(x, f, h, point[; method=:three])

Find the derivative of a bookend point at either `:begin` or `:end` of dataset.
Acceptable values for `method` include {`:three`, `:five`}.
"""
function endpoint(x::Vector{T}, f::Vector{T}, h::T, point::Symbol;
        method::Symbol=:three)::Float64 where {T<:Float64}
    i = (point == :begin ? 1 : (point == :end ? length(x) : nothing))
    if isnothing(i)
        throw(ArgumentError("`point` must be either `:begin` or `:end`."))
    end
    @inbounds return if point == :begin
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

@inline function endpoint(x::AbstractVector{T}, f::AbstractVector{T}, h::T, point::Symbol;
        method::Symbol=:three) where {T<:Float64}
    return endpoint(float(collect(x)), float(collect(f)), h, point; method=method)
end

"""
    midpoint(x, f, h, point[; method=:three])

Find the derivative of some point within a dataset.
Acceptable values for `method` include {`:three`, `:five`, `:2nd`}.
"""
function midpoint(x::Vector{T}, f::Vector{T}, h::T, point::Int64;
        method::Symbol=:three)::Float64 where {T<:Float64}
    i = point
    @inbounds return if method == :three
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

@inline function midpoint(x::AbstractVector{T}, f::AbstractVector{T}, h::T, point::Int64;
        method::Symbol=:three) where {T<:Float64}
    return midpoint(float(collect(x)), float(collect(f)), h, point; method=method)
end

end
