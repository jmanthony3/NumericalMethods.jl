struct SquareMatrixError <: Exception
    dims
end

struct SameSizeError <: Exception
    l1::Int64
    l2::Int64
end

struct SystemOfEquationError <: Exception
    A::Matrix{Float64}
    x::Vector{Float64}
end

struct SymmetricError <: Exception
end

struct PositiveDefiniteError <: Exception
end

struct IntervalBoundsError <: Exception
    f::Function
    a::Float64
    b::Float64
end

Base.showerror(io::IO, e::SquareMatrixError)        = print(io, "Matrix must be square. Is ", e.dims)
Base.showerror(io::IO, e::SameSizeError)            = print(io, "Vectors must be same the length.", e.l1, " ≠ ", e.l2)
Base.showerror(io::IO, e::SystemOfEquationError)    = print(io, "Systems of equations are not of compatible shape. Matrix is ", size(e.A), " and vector is ", size(e.x))
Base.showerror(io::IO, e::SymmetricError)           = print(io, "Matrix must be symmetric.")
Base.showerror(io::IO, e::PositiveDefiniteError)    = print(io, "Matrix must be positive definite.")
Base.showerror(io::IO, e::IntervalBoundsError)      = print(io, "Interval bounds must yield opposite signs in function, f:= [f(a=$(e.a)) = $(e.f(e.a)), f(b=$(e.b)) = $(e.f(e.b))].")