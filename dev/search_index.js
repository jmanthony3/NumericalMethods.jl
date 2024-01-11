var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = NumericalMethods","category":"page"},{"location":"#NumericalMethods","page":"Home","title":"NumericalMethods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NumericalMethods.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Functions","page":"Home","title":"Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [NumericalMethods]","category":"page"},{"location":"#NumericalMethods.ODE","page":"Home","title":"NumericalMethods.ODE","text":"ODE(f, a, b, h, α, N)\n\nStructure of the boundary conditions to differential equation where f is the time derivative of the function to approximate.\n\nNotes\n\nMake sure the independent variable is the first argument of f!\n\n\n\n\n\n","category":"type"},{"location":"#NumericalMethods.clamped-Union{Tuple{T}, Tuple{T, T, T}} where T<:(AbstractVector{T} where T)","page":"Home","title":"NumericalMethods.clamped","text":"clamped(x, f, fp)\n\nThe bookend polynomials will have the same slope entering and exiting the interval as the derivative at the respective endpoint.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalMethods.endpoint-Union{Tuple{T}, Tuple{T, T, Real, Symbol}} where T<:(AbstractVector{T} where T)","page":"Home","title":"NumericalMethods.endpoint","text":"endpoint(x, f, h, point[, point_type=\"three\"])\n\nFind the derivative of a bookend point at either :begin or :end of dataset. Acceptable values for method include {:three, :five}.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalMethods.integrate-Tuple{Union{Function, AbstractVector{T} where T}, AbstractVector{T} where T}","page":"Home","title":"NumericalMethods.integrate","text":"integrate(f[; rule=:trapezoidal, tol=10^-3])\n\nNotes\n\nFind the definite integral by some composite numeric quadrature. f may be a function or range. The domain may be defined with a vector, x or on the interval [a, b] either by number of sub-intervals, n or step-size, h. rule accepts {:trapezoidal (default), :midpoint, :simpson13, :simpson38, :simpsonN}. Dataset may contain unevenly spaces points.\n\nReferences\n\nhttps://en.wikipedia.org/wiki/Simpson%27s_rule\n\n\n\n\n\n","category":"method"},{"location":"#NumericalMethods.ivp-Tuple{ODE}","page":"Home","title":"NumericalMethods.ivp","text":"ivp(obj::ODE[, tol=10^-3; method=:forward_euler])\n\nSolve obj according to method ∈ {:forward_euler (default), :backward_euler, :improved_euler, :modified_euler, :runge_kutta}.\n\nNotes\n\nEach method has an equivalent convenience function. E.g. ivp(obj; method=:runge_kutta) ≡ runge_kutta(obj).\n\n\n\n\n\n","category":"method"},{"location":"#NumericalMethods.lagrange-Union{Tuple{T}, Tuple{T, T}, Tuple{T, T, Union{Nothing, Integer}}} where T<:(AbstractVector{T} where T)","page":"Home","title":"NumericalMethods.lagrange","text":"lagrange()\n\nGiven a domain and range, construct a Lagrangian polynomial. Polynomial will quickly begin to oscillate for larger datasets.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalMethods.linearleastsquares-Union{Tuple{T}, Tuple{T, T, Integer}} where T<:(AbstractVector{T} where T)","page":"Home","title":"NumericalMethods.linearleastsquares","text":"linearleastsquares(x, f, n::Integer)\n\nConstruct a polynomial of some degree, n while minimizing the least squares error.\n\nNotes\n\nLeast squares error := E = sum_i=1^my_i - P_n(x_i)^2\n\nConstructed polynomial of the form: P(x) = a_nx^n + a_n - 1x^n - 1 + dots + a_1x + a_0\n\n\n\n\n\n","category":"method"},{"location":"#NumericalMethods.linearleastsquares-Union{Tuple{T}, Tuple{T, T, Symbol}} where T<:(AbstractVector{T} where T)","page":"Home","title":"NumericalMethods.linearleastsquares","text":"linearleastsquares(domain, f, type::Symbol)\n\nGiven a domain and range, yield the coefficients for an equation and the equation of the form y = ax^b.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalMethods.midpoint-Union{Tuple{T}, Tuple{T, T, Real, Integer}} where T<:(AbstractVector{T} where T)","page":"Home","title":"NumericalMethods.midpoint","text":"midpoint(x, f, h, point[, point_type=\"three\"])\n\nFind the derivative of some point within a dataset. Acceptable values for method include {:three, :five, :2nd}.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalMethods.natural-Union{Tuple{T}, Tuple{T, T}} where T<:(AbstractVector{T} where T)","page":"Home","title":"NumericalMethods.natural","text":"natural(x, f)\n\nThe bookend polynomials do not assume the slope entering and exiting the interval as the derivative at the respective endpoint.\n\n\n\n\n\n","category":"method"},{"location":"#NumericalMethods.newtondifference-Union{Tuple{T}, Tuple{T, T, Real}} where T<:(AbstractVector{T} where T)","page":"Home","title":"NumericalMethods.newtondifference","text":"newtondifference(x, f, α[; dir::Symbol=:auto])\n\nGiven a domain and range, construct some polynomial by Newton's Divided Difference. 'forward' or 'backward' construction. Will be chosen automatically if not specified.\n\nNotes\n\nDirection will be chosen if not specified. Polynomials best made with even spacing in domain; although, this is not completely necessary.\n\n\n\n\n\n","category":"method"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
