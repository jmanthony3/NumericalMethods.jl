var documenterSearchIndex = {"docs":
[{"location":"numericaldifferentiation/#Numerical-Differentiation","page":"Numerical Differentiation","title":"Numerical Differentiation","text":"","category":"section"},{"location":"numericaldifferentiation/","page":"Numerical Differentiation","title":"Numerical Differentiation","text":"CurrentModule = NumericalMethods","category":"page"},{"location":"numericaldifferentiation/#Functions","page":"Numerical Differentiation","title":"Functions","text":"","category":"section"},{"location":"numericaldifferentiation/","page":"Numerical Differentiation","title":"Numerical Differentiation","text":"Modules = [NumericalMethods]\nOrder   = [:function]\nPages   = [\"numericaldifferentiation.jl\"]","category":"page"},{"location":"numericaldifferentiation/#NumericalMethods.endpoint-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, T, Symbol}} where T<:Float64","page":"Numerical Differentiation","title":"NumericalMethods.endpoint","text":"endpoint(x, f, h, point[; method=:three])\n\nFind the derivative of a bookend point at either :begin or :end of dataset. Acceptable values for method include {:three, :five}.\n\n\n\n\n\n","category":"method"},{"location":"numericaldifferentiation/#NumericalMethods.midpoint-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, T, Int64}} where T<:Float64","page":"Numerical Differentiation","title":"NumericalMethods.midpoint","text":"midpoint(x, f, h, point[; method=:three])\n\nFind the derivative of some point within a dataset. Acceptable values for method include {:three, :five, :2nd}.\n\n\n\n\n\n","category":"method"},{"location":"numericaldifferentiation/#NumericalMethods.n1derivative-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, Int64}} where T<:Float64","page":"Numerical Differentiation","title":"NumericalMethods.n1derivative","text":"n1derivative(x, f, j[; n=nothing])\n\nThe general (n + 1)-point formula to approximate f' at point j.\n\nNotes\n\nIf n = nothing, then entire dataset used to construct n'th Lagrange coefficient.\n\n\n\n\n\n","category":"method"},{"location":"numericaldifferentiation/#Index","page":"Numerical Differentiation","title":"Index","text":"","category":"section"},{"location":"numericaldifferentiation/","page":"Numerical Differentiation","title":"Numerical Differentiation","text":"Modules = [NumericalMethods]\nOrder   = [:type, :function]\nPages   = [\"numericaldifferentiation.md\"]","category":"page"},{"location":"initialvalueproblem/#Initial-Value-Problem","page":"Initial-Value Problem","title":"Initial-Value Problem","text":"","category":"section"},{"location":"initialvalueproblem/","page":"Initial-Value Problem","title":"Initial-Value Problem","text":"CurrentModule = NumericalMethods","category":"page"},{"location":"initialvalueproblem/#Types","page":"Initial-Value Problem","title":"Types","text":"","category":"section"},{"location":"initialvalueproblem/","page":"Initial-Value Problem","title":"Initial-Value Problem","text":"Modules = [NumericalMethods]\nOrder   = [:type]\nPages   = [\"initialvalueproblem.jl\"]","category":"page"},{"location":"initialvalueproblem/#NumericalMethods.InitialValueProblem","page":"Initial-Value Problem","title":"NumericalMethods.InitialValueProblem","text":"InitialValueProblem(f, a, b, h, α, N)\n\nStructure of the boundary conditions to Initial-Value Problem (IVP) differential equation.\n\nNotes\n\nMake sure the independent variable (e.g. time) is the first argument of f!\n\n\n\n\n\n","category":"type"},{"location":"initialvalueproblem/#Functions","page":"Initial-Value Problem","title":"Functions","text":"","category":"section"},{"location":"initialvalueproblem/","page":"Initial-Value Problem","title":"Initial-Value Problem","text":"Modules = [NumericalMethods]\nOrder   = [:function]\nPages   = [\"initialvalueproblem.jl\"]","category":"page"},{"location":"initialvalueproblem/#NumericalMethods.solve-Tuple{InitialValueProblem}","page":"Initial-Value Problem","title":"NumericalMethods.solve","text":"solve(ivp::InitialValueProblem[; method=:forward_euler, tol=10^-3])\n\nSolve ivp according to method ∈ {:forward_euler (default), :backward_euler, :improved_euler, :modified_euler, :runge_kutta}.\n\nNotes\n\nEach method has an equivalent convenience function. E.g. solve(ivp; method=:runge_kutta) ≡ runge_kutta(ivp).\n\n\n\n\n\n","category":"method"},{"location":"initialvalueproblem/#Index","page":"Initial-Value Problem","title":"Index","text":"","category":"section"},{"location":"initialvalueproblem/","page":"Initial-Value Problem","title":"Initial-Value Problem","text":"Modules = [NumericalMethods]\nOrder   = [:type, :function]\nPages   = [\"initialvalueproblem.md\"]","category":"page"},{"location":"singlevariableiteration/#Single-Variable-Iteration","page":"Single-Variable Iteration","title":"Single-Variable Iteration","text":"","category":"section"},{"location":"singlevariableiteration/","page":"Single-Variable Iteration","title":"Single-Variable Iteration","text":"CurrentModule = NumericalMethods","category":"page"},{"location":"singlevariableiteration/#Types","page":"Single-Variable Iteration","title":"Types","text":"","category":"section"},{"location":"singlevariableiteration/","page":"Single-Variable Iteration","title":"Single-Variable Iteration","text":"Modules = [NumericalMethods]\nOrder   = [:type]\nPages   = [\"singlevariableiteration.jl\"]","category":"page"},{"location":"singlevariableiteration/#NumericalMethods.SingleVariableIteration","page":"Single-Variable Iteration","title":"NumericalMethods.SingleVariableIteration","text":"SingleVariableIteration(f, a, b, n, tol)\n\nGiven f(p) such that p ∈ [a, b], find the root of a single-variable equation in so many iterations, n within tolerance, tol.\n\n\n\n\n\n","category":"type"},{"location":"singlevariableiteration/#Functions","page":"Single-Variable Iteration","title":"Functions","text":"","category":"section"},{"location":"singlevariableiteration/","page":"Single-Variable Iteration","title":"Single-Variable Iteration","text":"Modules = [NumericalMethods]\nOrder   = [:function]\nPages   = [\"singlevariableiteration.jl\"]","category":"page"},{"location":"singlevariableiteration/#NumericalMethods.bisection-Tuple{SingleVariableIteration}","page":"Single-Variable Iteration","title":"NumericalMethods.bisection","text":"bisection(svi::SingleVariableIteration)\n\nRoot-finding method: f(x) = 0. Search for solution by halving the bounds such that a and b initially yield opposite signs in function.\n\nNotes\n\nRelying on the Intermediate Value Theorem (IVT), this is a bracketed, root-finding method. This method is rather slow to converge but will always converge to a solution; therefore, is a good starter method.\n\n\n\n\n\n","category":"method"},{"location":"singlevariableiteration/#NumericalMethods.false_position-Tuple{SingleVariableIteration, Float64, Float64}","page":"Single-Variable Iteration","title":"NumericalMethods.false_position","text":"false_position(svi::SingleVariableIteration, p0, p1)\n\nAttempt root-finding method with initial guesses such that p0 and p1 in [a, b] yield opposite signs in function.\n\nUse function with lowest slope!\n\nNotes\n\nSimilar to, but slower to converge than, the secant_method() by including a test to ensure solution is root-bracketed.\n\nSee fixed_point() for theorem.\n\n\n\n\n\n","category":"method"},{"location":"singlevariableiteration/#NumericalMethods.fixed_point-Tuple{SingleVariableIteration, Float64}","page":"Single-Variable Iteration","title":"NumericalMethods.fixed_point","text":"fixed_point(svi::SingleVariableIteration, p0)\n\nAttempt root-finding method with initial guess, p0 in [a, b] by solving the equation g(p) = p via f(p) - p = 0.\n\nUse function with lowest slope!\n\nNot root-bracketed.\n\nNotes\n\nTheorem:\n\nExistence of a fixed-point:  If g ∈ C[a,b] and g(x) ∈ C[a, b] for all x ∈ [a, b], then function, g has a fixed point, p ∈ [a, b].\nUniqueness of a fixed point:  If g'(x) exists on [a, b] and a positive constant, k < 1 exist with {|g'(x)| ≤ k | x ∈ (a, b)}, then there is exactly one fixed-point, p ∈ [a, b].\n\nConverges by mathcalO(textlinear) if g'(p) ≠ 0, and mathcalO(textquadratic) if g'(p) = 0 and g''(p) < M, where M = g''(ξ) that is the error function.\n\n\n\n\n\n","category":"method"},{"location":"singlevariableiteration/#NumericalMethods.maximum_iterations","page":"Single-Variable Iteration","title":"NumericalMethods.maximum_iterations","text":"maximum_iterations(obj, method[, p0=0])\n\nFind greatest integer for maximum iterations within tolerance.\n\nNotes\n\nAcceptable values for method ∈ {:bisection, :fixed_point, :newton_raphson, :secant_method, :false_position}. Initial guess, p0 for function solution is not needed for :bisection method.\n\n\n\n\n\n","category":"function"},{"location":"singlevariableiteration/#NumericalMethods.maximum_slope-Tuple{SingleVariableIteration}","page":"Single-Variable Iteration","title":"NumericalMethods.maximum_slope","text":"maximum_slope(svi::SingleVariableIteration)\n\nFind the greatest value for first derivative of function.\n\n\n\n\n\n","category":"method"},{"location":"singlevariableiteration/#NumericalMethods.newton_raphson-Tuple{SingleVariableIteration, Float64}","page":"Single-Variable Iteration","title":"NumericalMethods.newton_raphson","text":"newton_raphson(svi::SingleVariableIteration, p0[; df=nothing])\n\nAttempt root-finding method with initial guess, p0 in [a, b] by solving the equation g(p) = p via f(p) - p = 0. df will be the first derivative of function if not given.\n\nUse function with lowest slope!\n\ndf(x) ≠ 0\n\nNotes\n\nQuickest convergence rate, but not root-bracketed and has trouble with symmetric functions! Initial guess, p0 must be close to real solution; else, will converge to different root or oscillate (if symmetric). This method can be viewed as fixed-point iteration.\n\nTechnique based on first Taylor polynomial expansion of f about p_0 (that is p0) and evaluated at x = p. |p - p_0| is assumed small; therefore, 2ⁿᵈ-order Taylor term, the error, is small.\n\nSee fixed_point() for theorem.\n\n\n\n\n\n","category":"method"},{"location":"singlevariableiteration/#NumericalMethods.secant_method-Tuple{SingleVariableIteration, Float64, Float64}","page":"Single-Variable Iteration","title":"NumericalMethods.secant_method","text":"secant_method(svi::SingleVariableIteration, p0, p1)\n\nAttempt root-finding method with initial guesses such that p0 and p1 in [a, b] yield opposite signs in function.\n\nUse function with lowest slope!\n\nNot root-bracketed.\n\nNotes\n\nMethod is less computationally expensive than newton_raphson() but may converge at slower rate by circumventing need to calculate derivative.\n\nSee fixed_point() for theorem.\n\n\n\n\n\n","category":"method"},{"location":"singlevariableiteration/#NumericalMethods.solve-Union{Tuple{SingleVariableIteration}, Tuple{T}} where T<:Float64","page":"Single-Variable Iteration","title":"NumericalMethods.solve","text":"solve(svi::SingleVariableIteration[; method=:bisection, p0, p1, df])\n\nAttempt to find where f(p) = 0 according to method ∈ {:bisection (default), :fixed_point, :newton_raphson, :secant_method, :false_position}. Each method also has a convenience function.\n\nNotes\n\nConvergence Rates:     :newton_raphson > :secant_method > :false_position > :fixed_point > :bisection\n\n:bisection is default because will always converge to a solution but is not the quickest.\n\n\n\n\n\n","category":"method"},{"location":"singlevariableiteration/#Index","page":"Single-Variable Iteration","title":"Index","text":"","category":"section"},{"location":"singlevariableiteration/","page":"Single-Variable Iteration","title":"Single-Variable Iteration","text":"Modules = [NumericalMethods]\nOrder   = [:type, :function]\nPages   = [\"singlevariableiteration.md\"]","category":"page"},{"location":"multivariableiteration/#Multi-Variable-Iteration","page":"Multi-Variable Iteration","title":"Multi-Variable Iteration","text":"","category":"section"},{"location":"multivariableiteration/","page":"Multi-Variable Iteration","title":"Multi-Variable Iteration","text":"CurrentModule = NumericalMethods","category":"page"},{"location":"multivariableiteration/#Types","page":"Multi-Variable Iteration","title":"Types","text":"","category":"section"},{"location":"multivariableiteration/","page":"Multi-Variable Iteration","title":"Multi-Variable Iteration","text":"Modules = [NumericalMethods]\nOrder   = [:type]\nPages   = [\"multivariableiteration.jl\"]","category":"page"},{"location":"multivariableiteration/#NumericalMethods.MultiVariableIteration","page":"Multi-Variable Iteration","title":"NumericalMethods.MultiVariableIteration","text":"MultiVariableIteration(A, x, b, N, tol)\n\nGiven the system of equations, mathbfAvecx = vecb, attempt to solve vecx = mathbfA^-1vecb in so many iterations, N within tolerance, tol. Ideal for large, sparse systems.\n\n\n\n\n\n","category":"type"},{"location":"multivariableiteration/#Functions","page":"Multi-Variable Iteration","title":"Functions","text":"","category":"section"},{"location":"multivariableiteration/","page":"Multi-Variable Iteration","title":"Multi-Variable Iteration","text":"Modules = [NumericalMethods]\nOrder   = [:function]\nPages   = [\"multivariableiteration.jl\"]","category":"page"},{"location":"multivariableiteration/#NumericalMethods.condition_number-Tuple{Matrix{T} where T}","page":"Multi-Variable Iteration","title":"NumericalMethods.condition_number","text":"condition_number(A)\n\nFind the condition number of matrix, A.\n\nNotes\n\nDefinition [burdenNumericalAnalysis2016]_: The condition number of the non-singular matrix, mathbfA relative to a norm, ||⋅|| is\n\n    K(mathbfA) = mathbfA  mathbfA^-1\n\nA matrix is well-conditioned if K(mathbfA) is close to 1 and is ill-conditioned if significantly greater than 1.\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#NumericalMethods.diagonality-Tuple{Matrix{T} where T}","page":"Multi-Variable Iteration","title":"NumericalMethods.diagonality","text":"diagonality(A)\n\nDetermines whether matrix, A is strictly, diagonally dominant.\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#NumericalMethods.find_omega","page":"Multi-Variable Iteration","title":"NumericalMethods.find_omega","text":"find_omega(mvi[, omega=0.])\n\nFind the optimum relaxation parameter, omega for mvi.A and mvi.x, if possible.\n\nNotes\n\nIf 0 < omega < 2, then method will converge regardless of choice for mvi.x. If matrix, mvi.A is tridiagonal, then the spectral radius of Gauss-Seidel's T-matrix, mathbfT_g is used to calculate ω = 2  (1 + sqrt(1 - spectral_radius(mathbfT_g)) [burdenNumericalAnalysis2016]_.\n\n\n\n\n\n","category":"function"},{"location":"multivariableiteration/#NumericalMethods.gauss_seidel-Tuple{MultiVariableIteration}","page":"Multi-Variable Iteration","title":"NumericalMethods.gauss_seidel","text":"gauss_seidel(mvi)\n\nSolve vecx = mathbfA^-1vecb via the Gauss-Seidel Method.\n\nNotes\n\nThis method improves on jacobi() by using the most recently calculated entries in the approximation vector, x at the end of each iteration. The core algorithm by which method marches through iterations:\n\n    vecx^(k) = bigl( (mathbfD - mathbfL)^-1 * mathbfU bigr) \n        vecx^(k - 1) + bigl( (mathbfD - mathbfL)^-1 bigr)  vecb\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#NumericalMethods.jacobi-Tuple{MultiVariableIteration}","page":"Multi-Variable Iteration","title":"NumericalMethods.jacobi","text":"jacobi(mvi)\n\nSolve vecx = mathbfA^-1vecb via the Jacobi Method to find vecx.\n\nNotes\n\nThe core algorithm by which method marches through iterations:\n\n    vecx^(k) = bigl( mathbfD^-1 * (mathbfL + mathbfU) bigr) \n        vecx^(k - 1) + ( mathbfD^-1 )  vecb\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#NumericalMethods.newton_raphson-Tuple{MultiVariableIteration, Tuple{Vararg{Symbolics.Num, N} where N}}","page":"Multi-Variable Iteration","title":"NumericalMethods.newton_raphson","text":"newton_raphson(mvi::MultiVariableIteration, variables::Tuple{Vararg{Num}})\n\nSolve non-linear systems of equations, vecx = mathbfA^-1vecb via the Newton-Raphson Method.\n\nHere, mvi.A should be a vector of functions wherein each variable is represented.\n\nExamples\n\nf1(x1, x2, x3)  = 3x1 - cos(x2*x3) - 0.5\nf2(x1, x2, x3)  = x1^2 - 81(x2 + 0.1)^2 + sin(x3) + 1.06\nf3(x1, x2, x3)  = exp(-x1*x2) + 20x3 + 3\\(10π - 3)\nA               = [f1, f2, f3]\nb               = zeros(length(A))\nx0              = [0.1, 0.1, -0.1]\ntol             = 1e-9\nmvi             = MultiVariableIteration(A, x0, b, 5, tol)\nusing Symbolics\n@variables x1, x2, x3\nnewton_raphson(mvi, (x1, x2, x3))\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#NumericalMethods.positive_definite-Tuple{Matrix{T} where T}","page":"Multi-Variable Iteration","title":"NumericalMethods.positive_definite","text":"positive_definite(A)\n\nDetermines whether matrix, A is positive definite.\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#NumericalMethods.solve-Tuple{MultiVariableIteration}","page":"Multi-Variable Iteration","title":"NumericalMethods.solve","text":"solve(mvi::MultiVariableIteration[; method=:jacobi, omega=0., variables])\n\nSolve vecx = mathbfA^-1vecb according to method ∈ {:jacobi (default), :gauss_seidel, :successive_relaxation, :newton_raphson}.\n\nEach method has an equivalent convenience function. E.g. solve(mvi; method=:jacobi) ≡ jacobi(ivp).\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#NumericalMethods.spectral_radius-Tuple{Matrix{T} where T}","page":"Multi-Variable Iteration","title":"NumericalMethods.spectral_radius","text":"spectral_radius(A)\n\nFinds the spectral radius of matrix, A.\n\nNotes\n\nρ(mathbfA) = maxλ, where λ is the set of eigvalsvalues for A [burdenNumericalAnalysis2016]_.\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#NumericalMethods.successive_relaxation","page":"Multi-Variable Iteration","title":"NumericalMethods.successive_relaxation","text":"successive_relaxation(mvi[, omega=0.])\n\nSolve vecx = mathbfA^-1vecb via the Successive Relaxation Method. Method is Successive Over-Relaxation (SOR) if omega > 1, Successive Under-Relaxation (SUR) if omega < 1, and reduces to Gauss-Seidel if omega = 1.\n\nNotes\n\nSOR and SUR accelerate or deccelerate convergence of gauss_seidel(), respectively, by decreasing or increasing the spectral radius of mvi.A. The core algorithm by which method marches through iterations:\n\n    vecx^(k) = bigl( (mathbfD - ωmathbfL)^-1 * ((1 - ω)*mathbfD + ωmathbfU) bigr) \n        vecx^(k - 1) + ω( (mathbfD - ωmathbfL)^-1 )  vecb\n\nIf left unspecified, and if possible, an optimum relaxation parameter, ω will be calculated by find_omega().\n\n\n\n\n\n","category":"function"},{"location":"multivariableiteration/#NumericalMethods.symmetry-Tuple{Matrix{T} where T}","page":"Multi-Variable Iteration","title":"NumericalMethods.symmetry","text":"symmetry(A)\n\nDetermines whether matrix, A is symmetric.\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#NumericalMethods.tridiagonality-Tuple{Matrix{T} where T}","page":"Multi-Variable Iteration","title":"NumericalMethods.tridiagonality","text":"tridiagonality(A)\n\nDetermine whether matrix, A is tridiagonal.\n\n\n\n\n\n","category":"method"},{"location":"multivariableiteration/#Index","page":"Multi-Variable Iteration","title":"Index","text":"","category":"section"},{"location":"multivariableiteration/","page":"Multi-Variable Iteration","title":"Multi-Variable Iteration","text":"Modules = [NumericalMethods]\nOrder   = [:type, :function]\nPages   = [\"multivariableiteration.md\"]","category":"page"},{"location":"numericalintegration/#Numerical-Integration","page":"Numerical Integration","title":"Numerical Integration","text":"","category":"section"},{"location":"numericalintegration/","page":"Numerical Integration","title":"Numerical Integration","text":"CurrentModule = NumericalMethods","category":"page"},{"location":"numericalintegration/#Functions","page":"Numerical Integration","title":"Functions","text":"","category":"section"},{"location":"numericalintegration/","page":"Numerical Integration","title":"Numerical Integration","text":"Modules = [NumericalMethods]\nOrder   = [:function]\nPages   = [\"numericalintegration.jl\"]","category":"page"},{"location":"numericalintegration/#NumericalMethods.integrate-Union{Tuple{T}, Tuple{Union{AbstractVector{T}, Function}, AbstractVector{T}}} where T<:Real","page":"Numerical Integration","title":"NumericalMethods.integrate","text":"integrate(f, x      [; rule=:trapezoidal, tol=10^-3])\nintegrate(f, a, b, h[; rule=:trapezoidal, tol=10^-3])\nintegrate(f, a, b, n[; rule=:trapezoidal, tol=10^-3])\n\nFind the definite integral by some numerical quadrature.\n\nNotes\n\nf may be a function or vector. The domain may be defined with a vector, x or on the interval [a, b] either by number of sub-intervals, n or step-size, h. rule accepts {:trapezoidal (default), :midpoint, :simpson13, :simpson38, :simpsonN}. Dataset may contain unevenly spaces points.\n\nReferences\n\nhttps://en.wikipedia.org/wiki/Simpson%27s_rule\n\n\n\n\n\n","category":"method"},{"location":"numericalintegration/#Index","page":"Numerical Integration","title":"Index","text":"","category":"section"},{"location":"numericalintegration/","page":"Numerical Integration","title":"Numerical Integration","text":"Modules = [NumericalMethods]\nOrder   = [:type, :function]\nPages   = [\"numericalintegration.md\"]","category":"page"},{"location":"interpolation/#Interpolation","page":"Interpolation","title":"Interpolation","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"CurrentModule = NumericalMethods","category":"page"},{"location":"interpolation/#Functions","page":"Interpolation","title":"Functions","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"Modules = [NumericalMethods]\nOrder   = [:function]\nPages   = [\"interpolation.jl\"]","category":"page"},{"location":"interpolation/#NumericalMethods.bezier-Union{Tuple{T}, NTuple{4, T}} where T<:(AbstractVector{T} where T)","page":"Interpolation","title":"NumericalMethods.bezier","text":"bezier(x, y, xguides, yguides)\n\nAn application of Hermitic polynomials to draw Bezier curves between points.\n\nNotes\n\nEach argument should be a one-to-one mapping of points, (xᵢ, yᵢ) and (xᵢ₊₁, yᵢ₊₁) and their respective guide points, (xᵢ⁺, yᵢ⁺) and (xᵢ₊₁⁻, yᵢ₊₁⁻).\n\n\n\n\n\n","category":"method"},{"location":"interpolation/#NumericalMethods.clamped-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, Vector{T}}} where T<:Float64","page":"Interpolation","title":"NumericalMethods.clamped","text":"clamped(x, f, fp)\n\nThe bookend polynomials will have the same slope entering and exiting the interval as the derivative at the respective endpoint.\n\n\n\n\n\n","category":"method"},{"location":"interpolation/#NumericalMethods.lagrange-Union{Tuple{T}, Tuple{T, T}} where T<:Vector{Float64}","page":"Interpolation","title":"NumericalMethods.lagrange","text":"lagrange(x, f[; n=nothing])\n\nGiven a domain, x and range, f, construct the nth Lagrangian polynomial.\n\nNotes\n\nIf n=nothing, then method will utilize entire dataset. Polynomials will quickly oscillate for larger datasets.\n\n\n\n\n\n","category":"method"},{"location":"interpolation/#NumericalMethods.linearleastsquares-Union{Tuple{T}, Tuple{T, T, Int64}} where T<:(AbstractVector{T} where T)","page":"Interpolation","title":"NumericalMethods.linearleastsquares","text":"linearleastsquares(x, f, n::Int64)\n\nConstruct a polynomial of degree, n while minimizing the least squares error.\n\nNotes\n\nLeast squares error := E = sum_i=1^my_i - P_n(x_i)^2\n\nConstructed polynomial of the form: P(x) = a_nx^n + a_n - 1x^n - 1 + dots + a_1x + a_0\n\n\n\n\n\n","category":"method"},{"location":"interpolation/#NumericalMethods.linearleastsquares-Union{Tuple{T}, Tuple{T, T, Symbol}} where T<:(AbstractVector{T} where T)","page":"Interpolation","title":"NumericalMethods.linearleastsquares","text":"linearleastsquares(x, f, type::Symbol)\n\nGiven a domain and range, yield the coefficients for an equation and the equation of the form y = ax^b.\n\n\n\n\n\n","category":"method"},{"location":"interpolation/#NumericalMethods.natural-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}}} where T<:Float64","page":"Interpolation","title":"NumericalMethods.natural","text":"natural(x, f)\n\nThe bookend polynomials do not assume the slope entering and exiting the interval as the derivative at the respective endpoint.\n\n\n\n\n\n","category":"method"},{"location":"interpolation/#NumericalMethods.newtondifference-Union{Tuple{T}, Tuple{Vector{T}, Vector{T}, Float64}} where T<:Float64","page":"Interpolation","title":"NumericalMethods.newtondifference","text":"newtondifference(x, f, α[; dir::Symbol=:auto])\n\nGiven a domain, x and range, f, construct some polynomial by Newton's Divided Difference centered around α. :forward or :backward construction.\n\nNotes\n\nDirection will be chosen if not specified. Polynomials best made with even spacing in x; although, this is not completely necessary.\n\n\n\n\n\n","category":"method"},{"location":"interpolation/#Index","page":"Interpolation","title":"Index","text":"","category":"section"},{"location":"interpolation/","page":"Interpolation","title":"Interpolation","text":"Modules = [NumericalMethods]\nOrder   = [:type, :function]\nPages   = [\"interpolation.md\"]","category":"page"},{"location":"#NumericalMethods","page":"Home","title":"NumericalMethods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for NumericalMethods.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NumericalMethods]\nPages   = [\n    \"interpolation.md\",\n    \"singlevariableiteration.md\",\n    \"numericaldifferentiation.md\",\n    \"numericalintegration.md\",\n    \"initialvalueproblem.md\",\n    \"multivariableiteration.md\"\n]\nDepth   = 1","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [NumericalMethods]\nPages   = [\n    \"interpolation.md\",\n    \"singlevariableiteration.md\",\n    \"numericaldifferentiation.md\",\n    \"numericalintegration.md\",\n    \"initialvalueproblem.md\",\n    \"multivariableiteration.md\"\n]","category":"page"}]
}
