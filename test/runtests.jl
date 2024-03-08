using NumericalMethods
using Test

@testset "NumericalMethods.jl" begin
    # Write your tests here.
    ## single-variable iteration
    f(x)    = x^3 + 4x^2 - 10
    a, b, N, tol = 1., 2., 50, 10^-9
    svi     = SingleVariableIteration(f, a, b, N, tol)
    @test round(bisection(svi),                     digits=6)       == 1.365230
    g(x)    = 2 \ √(10 - x^3)
    svi     = SingleVariableIteration(g, a, b, N, tol)
    @test round(fixed_point(svi, 1.5),              digits=6)       == 1.365230
    h(x)    = cos(x) - x
    a, b, N, tol = 0., pi/2., 50, 10^-6
    svi     = SingleVariableIteration(h, a, b, N, tol)
    @test round(newton_raphson(svi, pi / 4.),       digits=6)       == 0.739085
    @test round(secant_method(svi, 0.5, pi / 4.),   digits=6)       == 0.739085
    a, b, N, tol = 0., pi/2., 100, 10^-1
    svi     = SingleVariableIteration(h, a, b, N, tol)
    @test round(false_position(svi, 0.5, pi / 4.),  digits=6)       == 0.739058

    ## interpolation
    x       = [0.01, 0.15, 0.31, 0.5, 0.6, 0.75]
    y       = [1.0, 1.004, 1.031, 1.117, 1.223, 1.422]
    n       = 2
    @test round(linearleastsquares(x, y, n)[2],     digits=6)       == 0.000873

    a, b, h = 1., 2.2, 0.3
    x       = [a:h:b...]
    y       = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
    @test round(newtondifference(x, y, 1.5)(1.5),   digits=6)       == 0.511820

    x       = float.([0,  5,  9, 12, 16, 23, 28])   # days
    y1      = float.([5, 14, 39, 34, 28, 26, 25])   # mg
    y2      = float.([5, 13, 15, 14, 12, 11, 10])   # mg
    @test round(lagrange(x, y1)[2][end],            digits=6)       == 0.006033
    @test round(lagrange(x, y2)[2][end],            digits=6)       == 0.000283

    x       = [0.0:0.1:0.5...]
    y       = sin.(3 .* x)
    dx      = 3cos.(3 .* x)
    @test round(clamped(x, y, dx)[2][begin](0.05),  digits=6)       == 0.149437
    @test round(natural(x, y)[2][begin](0.05),      digits=6)       == 0.149408

    ## derivatives
    h, a, b = 0.25, 1, 2
    x       = [a:h:b...]
    y       = (x .^ 2) .* cos.(x)
    @test round(n1derivative(x, y, 2),              digits=6)       == -0.690759
    @test round(endpoint(x, y, h, :begin),          digits=6)       ==  0.381398
    @test round(midpoint(x, y, h, 2),               digits=6)       == -0.762287

    ## integration
    @test round(integrate(y, x, rule=:trapezoidal), digits=6)       == -0.114043
    @test round(integrate(y, x, rule=:simpson13),   digits=6)       == -0.084893
    @test round(integrate(y, x, rule=:simpson38),   digits=6)       == -0.213298
    @test round(integrate(y, x, rule=:simpsonN),    digits=6)       == -0.084893
    @test round(integrate(y, x, rule=:midpoint),    digits=6)       == -0.017729

    ## ODE/PDE
    f(t, y) = y - t^2 + 1
    a, h    = 0., 0.2
    N, α    = 10, 0.5
    ivp     = InitialValueProblem(f, a, h, α, N)
    @test round.(runge_kutta(ivp),                  digits=6)       == [
        0.5,        0.829293,   1.214076,
        1.648922,   2.127203,   2.640823,
        3.179894,   3.732340,   4.283409,
        4.815086,   5.305363]
end
