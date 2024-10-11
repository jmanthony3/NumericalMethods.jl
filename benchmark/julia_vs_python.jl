using BenchmarkTools
using NumericalMethods
using PyCall
using Suppressor

nm = pyimport("joby_m_anthony_iii.numerical_methods")

A   = [10. -1. 2. 0.; -1. 11. -1. 3.; 2. -1. 10. -1.; 0. 3. -1. 8.]
b   = [6., 25., -11., 15]
x0  = [0., 0., 0., 0.]
tol = 1e-3

obj = @suppress nm.MultiVariableIteration(A, x0, b, -3, 10)
println("Jacobi (Python): ", @btime @suppress begin
    ($obj).jacobi()["Approximations"][end]
end)

MVI = MultiVariableIteration(A, x0, b, 10, tol)
println("Jacobi (Julia): ", @btime jacobi($MVI))