module NumericalMethods

include("base.jl")
export solve
export newton_raphson
export diagonality
export spectral_radius
export condition_number
export symmetry
export positive_definite
export tridiagonality

include("singlevariableiteration.jl")
using .SingleVariableIterations
export SingleVariableIteration
# export solve
export bisection
export fixed_point
# export newton_raphson
export secant_method
export false_position

include("interpolation.jl")
using .Interpolations
export linearinterpolation
export lagrange
export newtondifference
export natural
export clamped
export bezier
export linearleastsquares

include("numericaldifferentiation.jl")
using .Derivatives
export n1derivative
export endpoint
export midpoint

include("numericalintegration.jl")
using .Integrations
export integrate

include("initialvalueproblem.jl")
using .InitialValueProblems
export InitialValueProblem
# export solve
export forward_euler
export backward_euler
export improved_euler
export modified_euler
export runge_kutta

include("multivariableiteration.jl")
using .MultiVariableIterations
export MultiVariableIteration
# export diagonality
# export spectral_radius
# export condition_number
# export symmetry
# export positive_definite
# export tridiagonality
export find_omega
export jacobian_form
# export solve
export gauss_seidel
export jacobi
# export newton_raphson
export successive_relaxation

include("systemofequation.jl")
using .SystemOfEquations
export SystemOfEquation
export conjugate_gradient
export gaussian_elimination
export steepest_descent
# export solve

include("eigenvalue.jl")
using .Eigenvalues
export Eigenvalue
export power_method
export inverse_power_method
export qr_algorithm

end
