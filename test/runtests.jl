using LinearAlgebra: diag, diagm
using NumericalMethods
using Symbolics: @variables # KEY: [10.1145/3511528.3511535]
using Test

@testset verbose=true "NumericalMethods.jl" begin
    @testset "Single-Variable Iteration" begin
        # single-variable iteration
        f(x)    = x^3 + 4x^2 - 10
        a, b, N, tol = 1., 2., 50, 10^-9
        SVI     = SingleVariableIteration(f, a, b, N, tol)
        @test round(bisection(SVI);                     digits=6)       == 1.365230
        g(x)    = 2 \ √(10 - x^3)
        SVI     = SingleVariableIteration(g, a, b, N, tol)
        @test round(fixed_point(SVI, 1.5);              digits=6)       == 1.365230
        h(x)    = cos(x) - x
        a, b, N, tol = 0., pi/2., 50, 10^-6
        SVI     = SingleVariableIteration(h, a, b, N, tol)
        @test round(newton_raphson(SVI, pi / 4.);       digits=6)       == 0.739085
        @test round(secant_method(SVI, 0.5, pi / 4.);   digits=6)       == 0.739085
        a, b, N, tol = 0., pi/2., 100, 10^-1
        SVI     = SingleVariableIteration(h, a, b, N, tol)
        @test round(false_position(SVI, 0.5, pi / 4.);  digits=6)       == 0.739058
    end
    @testset "Interpolation" begin
        # interpolation
        x       = [0.01, 0.15, 0.31, 0.5, 0.6, 0.75]
        y       = [1.0, 1.004, 1.031, 1.117, 1.223, 1.422]
        n       = 2
        @test round(linearleastsquares(x, y, n)[2];     digits=6)       == 0.000873

        a, b, h = 1., 2.2, 0.3
        x       = [a:h:b...]
        y       = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
        @test round(newtondifference(x, y, 1.5)(1.5);   digits=6)       == 0.511820

        x       = float.([0,  5,  9, 12, 16, 23, 28])   # days
        y1      = float.([5, 14, 39, 34, 28, 26, 25])   # mg
        y2      = float.([5, 13, 15, 14, 12, 11, 10])   # mg
        @test round(lagrange(x, y1)[2][end];            digits=6)       == 0.006033
        @test round(lagrange(x, y2)[2][end];            digits=6)       == 0.000283

        x       = [0.0:0.1:0.5...]
        y       = sin.(3 .* x)
        dx      = 3cos.(3 .* x)
        @test round(clamped(x, y, dx)[2][begin](0.05);  digits=6)       == 0.149437
        @test round(natural(x, y)[2][begin](0.05);      digits=6)       == 0.149408
    end
    @testset "Numerical Differentiation" begin
        # derivatives
        h, a, b = 0.25, 1, 2
        x       = [a:h:b...]
        y       = (x .^ 2) .* cos.(x)
        @test round(n1derivative(x, y, 2);              digits=6)       == -0.690759
        @test round(endpoint(x, y, h, :begin);          digits=6)       ==  0.381398
        @test round(midpoint(x, y, h, 2);               digits=6)       == -0.762287
    end
    @testset "Numerical Integration" begin
        # integration
        h, a, b = 0.25, 1, 2
        x       = [a:h:b...]
        y       = (x .^ 2) .* cos.(x)
        @test round(integrate(y, x, rule=:trapezoidal); digits=6)       == -0.114043
        @test round(integrate(y, x, rule=:simpson13);   digits=6)       == -0.084893
        @test round(integrate(y, x, rule=:simpson38);   digits=6)       == -0.213298
        @test round(integrate(y, x, rule=:simpsonN);    digits=6)       == -0.084893
        @test round(integrate(y, x, rule=:midpoint);    digits=6)       == -0.017729
    end
    @testset "Initial-Value Problem" begin
        # IVP
        f(t, y) = y - t^2 + 1
        a, h    = 0., 0.2
        N, α    = 10, 0.5
        IVP     = InitialValueProblem(f, a, h, α, N)
        @test round.(runge_kutta(IVP);                  digits=6)       == [
            0.5,        0.829293,   1.214076,
            1.648922,   2.127203,   2.640823,
            3.179894,   3.732340,   4.283409,
            4.815086,   5.305363]
    end
    @testset "Multi-Variable Iteration" begin
        # multi-variable iteration
        ## jacobi ~ gauss_seidel
        A   = [10. -1. 2. 0.; -1. 11. -1. 3.; 2. -1. 10. -1.; 0. 3. -1. 8.]
        b   = [6., 25., -11., 15]
        x0  = [0., 0., 0., 0.]
        tol = 1e-3
        MVI = MultiVariableIteration(A, x0, b, 10, tol)
        @test round.(jacobi(MVI);                       digits=4)       == [
            1.0001,     1.9998,    -0.9998,     0.9998]
        MVI = MultiVariableIteration(A, x0, b, 5, tol)
        @test round.(gauss_seidel(MVI);                 digits=4)       == [
            1.0001,     2.0000,    -1.0000,     1.0000]
        ## gauss_seidel ~ successive_relaxation
        A   = [4. 3. 0.; 3. 4. -1.; 0. -1. 4.]
        b   = [24., 30., -24.]
        x0  = ones(length(b))
        tol = 1e-3
        omega = 1.25
        MVI = MultiVariableIteration(A, x0, b, 7, tol)
        @test round.(gauss_seidel(MVI);                 digits=7)       == [
            3.0134110,  3.9888241, -5.0027940]
        @test all(isapprox.(round.(successive_relaxation(MVI, omega); digits=7), [
            3.0000498,  4.0002586, -5.0003486]; atol=10tol))
        ## newton_raphson
        f1(x1, x2, x3) = 3x1 - cos(x2*x3) - 0.5
        f2(x1, x2, x3) = x1^2 - 81(x2 + 0.1)^2 + sin(x3) + 1.06
        f3(x1, x2, x3) = exp(-x1*x2) + 20x3 + 3\(10π - 3)
        @variables x1, x2, x3
        variables = (x1, x2, x3)
        A   = [f1, f2, f3]
        b   = zeros(length(A))
        x0  = [0.1, 0.1, -0.1]
        tol = 1e-9
        MVI = MultiVariableIteration(A, x0, b, 5, tol)
        @test all(isapprox.(round.(newton_raphson(MVI, variables); digits=10), [
            0.5000000000, -1.375e-11, -0.5235987756]; atol=10tol))
    end
    @testset "System of Equations" begin
        # systems of equations
        ## gaussian_elimination
        A = float.([
            [1 1 0 3]
            [2 1 -1 1]
            [3 -1 -1 2]
            [-1 2 3 -1]
        ])
        b = float.([4, 1, -3, 4])
        tol = 10^-3.
        SOE = SystemOfEquation(A, b, 100, tol)
        @test all(isapprox.(round.(gaussian_elimination(SOE);                           digits=1), [
            -1, 2, 0, 1]; atol=10tol))
        ## steepest_descent ~ gaussian_elimination
        A   = float.([
            [4 1 1 0 1]
            [1 3 1 1 0]
            [1 1 5 -1 -1]
            [0 1 -1 4 0]
            [1 0 -1 0 4]])
        b   = fill(6., size(A)[1])
        tol = 10^-5.
        x0  = ones(length(b))
        SOE = SystemOfEquation(A, b, 100, tol)
        @test all(isapprox.(round.(steepest_descent(SOE, x0);                           digits=6), [
            0.451611, 0.709681, 1.677415, 1.741933, 1.806451]; atol=10tol))
        @test all(isapprox.(round.(conjugate_gradient(SOE, x0);                         digits=6), [
            0.451611, 0.709681, 1.677415, 1.741933, 1.806451]; atol=10tol))
        @test all(isapprox.(round.(conjugate_gradient(SOE, x0, diagm(0 => diag(A)));    digits=6), [
            0.451611, 0.709681, 1.677415, 1.741933, 1.806451]; atol=10tol))
    end
    @testset "Eigenvalue" begin
        # eigenvalue
        ## (inverse_)power_method
        A = float.([
            [-4 14 0]
            [-5 13 0]
            [-1 0 2]])
        tol = 10^-4
        x0 = ones(size(A)[1])
        EV = Eigenvalue(A, 100, tol)
        @test isapprox(round(power_method(EV, x0);                  digits=6), 6.000837; atol=10tol)
        @test isapprox(round(inverse_power_method(EV, x0, 19/3);    digits=6), 6.000017; atol=10tol)

        ## qr_algorithm
        A = float.([
            [3 1 0]
            [1 3 1]
            [0 1 3]])
        tol = 10^-5
        @test all(isapprox.(round.(qr_algorithm(Eigenvalue(A, 100, tol)); digits=5), [
            4.41420, 3.00000, 1.58579]; atol=10tol))
    end
    @testset "Boundary-Value Problem" begin
        # BVP
        ## linear_shooting_method ~ finite_difference_method
        p(x) = -x\2
        q(x) = 2/(x^2)
        r(x) = sin(log(x))/x^2
        BVP = BoundaryValueProblem([p, q, r], 1., 2., 0.1, 1., 2., 10)
        tol = 10^-6
        @test all(isapprox.(round.([y[1] for y in linear_shooting_method(BVP)]; digits=6), [
            1.0, 1.092629, 1.187085, 1.283382, 1.381446, 1.481159, 1.582392, 1.685014,
            1.788898, 1.893930, 2.0]; atol=10tol))
        @test all(isapprox.(round.(finite_difference_method(BVP; tol=tol); digits=6), [
            1.0, 1.092600, 1.187043, 1.283337, 1.381402, 1.481120, 1.582360, 1.684989,
            1.788882, 1.893921, 2.0]; atol=10tol))
    end
end