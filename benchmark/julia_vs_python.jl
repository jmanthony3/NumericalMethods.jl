using BenchmarkTools
using NumericalMethods
using PyCall
using Suppressor

nm = pyimport("joby_m_anthony_iii.numerical_methods")

py"""
import logging

# store the current log level to restore it later
original_log_level = logging.getLogger().getEffectiveLevel()
# set the log level to a higher level, e.g., WARNING or CRITICAL
logging.disable(logging.CRITICAL)
"""

println("# Comparing Python and Julia Implementations")

# Example 1 (p.456) of burden2015numerical
println("\nExample 1 on p. 456 of Burden 2015...")
A   = [10. -1. 2. 0.; -1. 11. -1. 3.; 2. -1. 10. -1.; 0. 3. -1. 8.]
b   = [6., 25., -11., 15]
x0  = [0., 0., 0., 0.]
tol = 1e-3

obj = nm.MultiVariableIteration(A, x0, b, -3, 10)
println("> Jacobi (Python): ")
@btime ($obj).jacobi()["Approximations"][end];

MVI = MultiVariableIteration(A, x0, b, 10, tol)
println("> Jacobi (Julia): ")
@btime jacobi($MVI);


println("\nOn matrices of increasing size with random elements...")
for n in (10, 25, 50, 100)
    local A   = rand(Float64, n, n)
    local b   = rand(Float64, n)
    local x0  = zeros(n)
    local N   = 10
    local tol = 1e-3

    local obj = nm.MultiVariableIteration(A, x0, b, -3, N)
    println("> Jacobi (Python) [$(n)x$(n)]: ")
    @btime ($obj).jacobi()["Approximations"][end];

    local MVI = MultiVariableIteration(A, x0, b, N, tol)
    println("> Jacobi (Julia) [$(n)x$(n)]: ")
    @btime jacobi($MVI);
end

# restore the original log level after the tests
py"logging.disable(original_log_level)"