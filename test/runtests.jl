using NumericalMethods
using StaticArrays
using Test

@testset "NumericalMethods.jl" begin
    # Write your tests here.
    ## derivatives
    h, a, b = 0.25, 1, 2
    x       = [a:h:b...]
    y       = (x .^ 2) .* cos.(x)
    n       = length(x)
    xstatic = SVector{n, AbstractFloat}(x)
    ystatic = SVector{n, AbstractFloat}(y)
    @test round(n1derivative(x, y, 2), digits=6)                    == -0.690759 # ... 4120392901
    @test round(n1derivative(xstatic, ystatic, 2), digits=6)        == -0.690759 # ... 4120392901
    @test round(endpoint(x, y, h, :begin), digits=6)                ==  0.381398 # ... 2872273567
    @test round(endpoint(xstatic, ystatic, h, :begin), digits=6)    ==  0.381398 # ... 2872273567
    @test round(midpoint(x, y, h, 2), digits=6)                     == -0.762287 # ... 2042316164
    @test round(midpoint(xstatic, ystatic, h, 2), digits=6)         == -0.762287 # ... 2042316164

    ## integration
    @test round(integrate(y, x, rule=:trapezoidal), digits=6)       == -0.114043 # ... 79264796141
    @test round(integrate(y, x, rule=:simpson13),   digits=6)       == -0.084893 # ... 08746263464
    @test round(integrate(y, x, rule=:simpson38),   digits=6)       == -0.213298 # ... 30449138726
    @test round(integrate(y, x, rule=:simpsonN),    digits=6)       == -0.084893 # ... 08746263458
    @test round(integrate(y, x, rule=:midpoint),    digits=6)       == -0.017729 # ... 11806132075

    ## interpolation
    x       = [0.01, 0.15, 0.31, 0.5, 0.6, 0.75]
    y       = [1.0, 1.004, 1.031, 1.117, 1.223, 1.422]
    degree  = 2
    @test round(linearleastsquares(x, y, degree)[2], digits=6)      == 0.000873 # ... 6264091126483

    a, b, h = 1., 2.2, 0.3
    x       = [a:h:b...]
    y       = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
    @test round(newtondifference(x, y, 1.5)(1.5), digits=6)         == 0.511820 # ... 9942386833
    n = length(x);
    x = SVector{n, AbstractFloat}(x);
    y = SVector{n, AbstractFloat}(y);
    @test round(newtondifference(x, y, 1.5)(1.5), digits=6)         == 0.511820 # ... 9942386833

    x       = float.([0,  5,  9, 12, 16, 23, 28])   # days
    y1      = float.([5, 14, 39, 34, 28, 26, 25])   # mg
    y2      = float.([5, 13, 15, 14, 12, 11, 10])   # mg
    @test lagrange(x, y1)[2][end]                                   == 0.
    @test lagrange(x, y2)[2][end]                                   == 0.
    n       = length(x)
    xstatic = SVector{n, AbstractFloat}(x)
    y1stat  = SVector{n, AbstractFloat}(y1)
    y2stat  = SVector{n, AbstractFloat}(y2)
    @test lagrange(xstatic, y1stat)[2][end]                         == 0.
    @test lagrange(xstatic, y2stat)[2][end]                         == 0.

    x   = [0.0:0.1:0.5...]
    y   = sin.(3 .* x)
    dx  = 3cos.(3 .* x)
    @test round(clamped(x, y, dx)[2][begin](0.05),  digits=6)       == 0.149437 # ... 07092796263
    @test round(natural(x, y)[2][begin](0.05),      digits=6)       == 0.149408 # ... 86040416965
end
