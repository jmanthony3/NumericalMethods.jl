using BenchmarkTools
using LUSE_ENGR701_704_NumericalMethods
using Symbolics: @variables

SUITE = BenchmarkGroup()
SUITE["svi"] = BenchmarkGroup(["iteration"])
SUITE["interp"] = BenchmarkGroup(["interpolation"])
SUITE["diff"] = BenchmarkGroup(["quadrature"])
SUITE["int"] = BenchmarkGroup(["quadrature"])
SUITE["ode"] = BenchmarkGroup(["iteration", "quadrature", "soe"])
SUITE["mvi"] = BenchmarkGroup(["iteration", "soe"])

# single-variable iteration
    f_svi(x)= x^3 + 4x^2 - 10
    a, b, N, tol = 1., 2., 50, 10^-9
    svi     = SingleVariableIteration(f_svi, a, b, N, tol)
    SUITE["svi"]["bisection"]               = @benchmarkable bisection(svi);
    g_svi(x)= 2 \ √(10 - x^3)
    svi     = SingleVariableIteration(g_svi, a, b, N, tol)
    SUITE["svi"]["fixed_point"]             = @benchmarkable fixed_point(svi, 1.5);
    h_svi(x)= cos(x) - x
    a, b, N, tol = 0., pi/2., 50, 10^-6
    svi     = SingleVariableIteration(h_svi, a, b, N, tol)
    SUITE["svi"]["newton_raphson"]          = @benchmarkable newton_raphson(svi, pi / 4.);
    SUITE["svi"]["secant_method"]           = @benchmarkable secant_method(svi, 0.5, pi / 4.);
    a, b, N, tol = 0., pi/2., 100, 10^-1
    svi     = SingleVariableIteration(h_svi, a, b, N, tol)
    SUITE["svi"]["false_position"]          = @benchmarkable false_position(svi, 0.5, pi / 4.);



# interpolation
    x       = [0.01, 0.15, 0.31, 0.5, 0.6, 0.75]
    y       = [1.0, 1.004, 1.031, 1.117, 1.223, 1.422]
    n       = 2
    SUITE["interp"]["linearleastsquares"]   = @benchmarkable linearleastsquares(x, y, n)[2];

    a, b, h = 1., 2.2, 0.3
    x       = [a:h:b...]
    y       = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
    SUITE["interp"]["newtondifference"]     = @benchmarkable newtondifference(x, y, 1.5)(1.5);

    x       = float.([0,  5,  9, 12, 16, 23, 28])   # days
    y1      = float.([5, 14, 39, 34, 28, 26, 25])   # mg
    y2      = float.([5, 13, 15, 14, 12, 11, 10])   # mg
    SUITE["interp"]["lagrange"]             = @benchmarkable lagrange(x, y1)[2][end];
    # SUITE["interp"]["lagrange"]             = @benchmarkable lagrange(x, y2)[2][end];

    x       = [0.0:0.1:0.5...]
    y       = sin.(3 .* x)
    dx      = 3cos.(3 .* x)
    SUITE["interp"]["clamped"]              = @benchmarkable clamped(x, y, dx)[2][begin](0.05);
    SUITE["interp"]["natural"]              = @benchmarkable natural(x, y)[2][begin](0.05);



# derivatives
    h, a, b = 0.25, 1, 2
    x       = [a:h:b...]
    y       = (x .^ 2) .* cos.(x)
    SUITE["diff"]["n1derivative"]           = @benchmarkable n1derivative(x, y, 2);
    SUITE["diff"]["endpoint"]               = @benchmarkable endpoint(x, y, h, :begin);
    SUITE["diff"]["midpoint"]               = @benchmarkable midpoint(x, y, h, 2);



# integration
    h, a, b = 0.25, 1, 2
    x       = [a:h:b...]
    y       = (x .^ 2) .* cos.(x)
    SUITE["int"]["trapezoidal"]             = @benchmarkable integrate(y, x, rule=:trapezoidal);
    SUITE["int"]["simpson13"]               = @benchmarkable integrate(y, x, rule=:simpson13);
    SUITE["int"]["simpson38"]               = @benchmarkable integrate(y, x, rule=:simpson38);
    SUITE["int"]["simpsonN"]                = @benchmarkable integrate(y, x, rule=:simpsonN);
    SUITE["int"]["midpoint"]                = @benchmarkable integrate(y, x, rule=:midpoint);



# ODE/PDE
    f_ode(t, y) = y - t^2 + 1
    a, h    = 0., 0.2
    N, α    = 10, 0.5
    ivp     = InitialValueProblem(f_ode, a, h, α, N)
    SUITE["ode"]["runge_kutta"]             = @benchmarkable runge_kutta(ivp);



# multi-variable iteration
    ## jacobi ~ gauss_seidel
        A   = [10. -1. 2. 0.; -1. 11. -1. 3.; 2. -1. 10. -1.; 0. 3. -1. 8.]
        b   = [6., 25., -11., 15]
        x0  = [0., 0., 0., 0.]
        tol = 1e-3
        mvi = MultiVariableIteration(A, x0, b, 10, tol)
        # SUITE["mvi"]["jacobi"]              = @benchmarkable jacobi(mvi);
        # mvi = MultiVariableIteration(A, x0, b, 5, tol)
        # SUITE["mvi"]["gauss_seidel"]        = @benchmarkable gauss_seidel(mvi);
    ## gauss_seidel ~ successive_relaxation
        A   = [4. 3. 0.; 3. 4. -1.; 0. -1. 4.]
        b   = [24., 30., -24.]
        x0  = ones(length(b))
        tol = 1e-3
        omega = 1.25
        mvi = MultiVariableIteration(A, x0, b, 7, tol)
        # SUITE["mvi"]["gauss_seidel"]        = @benchmarkable gauss_seidel(mvi);
        # SUITE["mvi"]["successive_relaxation"]= @benchmarkable successive_relaxation(mvi, omega);
    ## newton_raphson
        f_mvi1(x1, x2, x3) = 3x1 - cos(x2*x3) - 0.5
        f_mvi2(x1, x2, x3) = x1^2 - 81(x2 + 0.1)^2 + sin(x3) + 1.06
        f_mvi3(x1, x2, x3) = exp(-x1*x2) + 20x3 + 3\(10π - 3)
        @variables x1, x2, x3
        variables = (x1, x2, x3)
        A   = [f_mvi1, f_mvi2, f_mvi3]
        b   = zeros(length(A))
        x0  = [0.1, 0.1, -0.1]
        tol = 1e-9
        mvi = MultiVariableIteration(A, x0, b, 5, tol)
        # SUITE["mvi"]["newton_raphson"]      = @benchmarkable newton_raphson(mvi, variables);