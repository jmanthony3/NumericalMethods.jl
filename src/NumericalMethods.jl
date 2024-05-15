module NumericalMethods

export solve

export SingleVariableIteration
export bisection
export fixed_point
export newton_raphson
export secant_method
export false_position
include("singlevariableiteration.jl")

export linearinterpolation
export lagrange
export newtondifference
export natural
export clamped
export bezier
export linearleastsquares
include("interpolation.jl")

export n1derivative
export endpoint
export midpoint
include("numericaldifferentiation.jl")

export integrate
include("numericalintegration.jl")

export InitialValueProblem
export forward_euler
export backward_euler
export improved_euler
export modified_euler
export runge_kutta
include("initialvalueproblem.jl")

export MultiVariableIteration
export diagonality
export spectral_radius
export condition_number
export symmetry
export positive_definite
export tridiagonality
export find_omega
export jacobian_form
export gauss_seidel
export jacobi
export newton_raphson
export successive_relaxation
include("multivariableiteration.jl")

end
