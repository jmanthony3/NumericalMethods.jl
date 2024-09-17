# NumericalMethods

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jmanthony3.github.io/NumericalMethods.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jmanthony3.github.io/NumericalMethods.jl/dev/)
[![Build Status](https://github.com/jmanthony3/NumericalMethods.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jmanthony3/NumericalMethods.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Documentation](https://github.com/jmanthony3/NumericalMethods.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/jmanthony3/NumericalMethods.jl/actions/workflows/Documentation.yml)

This package provides:
- Native Julia implementation of some numerical methods from the 10ᵗʰ edition of _Numerical Analysis_ by Burden et al.
- `solve()` includes method dispatching for the following types:
  - `SingleVariableIteration`
  - `InitialValueProblem`
  - `MultiVariableIteration`

Each of the method dispatches on `solve()` include convenience functions for the respective numerical method as specified by the `method` keyword argument.
E.g. `solve(mvi::MultiVariableIteration; method=:jacobi) ≡ jacobi(mvi)` which applies the **Jacobi Iterative Technique** onto the system of equations defined in a `MultiVariableIteration` structure.

Key comments:
- Greatly inspired by the Python package [`joby_m_anthony_iii`](https://pypi.org/project/joby-m-anthony-iii/)
- Leverages [`StaticArrays.jl`](https://juliaarrays.github.io/StaticArrays.jl/stable/) (where appropriate) and practices of recommended [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)
- Heavy use of [`Symbolics.jl`]() functions--`build_function`, `simplify`, `expand_derivatives`, and `Differential`--for symbolic differentiation and evaluation:
  - `jacobian_form()`
  - `lagrange()`
  - `n1derivative()`
  - `maximum_slope()`
  - `newton_raphson()`
- `newton_raphson()` will solve quicker if the functional form of the derivative or Jacobian matrix is supplied

## Roadmap
- [x] Single-Variable Iteration (SVI)
  - [x] Bisection Method
  - [x] Fixed-Point Method
  - [x] Secant Method
  - [x] False Position
  - [x] Newton-Raphson
- [ ] Interpolation methods
  - [x] Cubic splines
  - [x] Newton's Divided Difference
  - [x] Lagrange Polynomials
  - [ ] Hermite Polynomials
  - [x] Linear Least Squares
  - [x] Bezier Curves
- [ ] Extrapolation methods
  - [ ] Richardson's Extrapolation
- [x] Solving for eigenvalues
  - [x] Power Method
  - [x] Inverse Power Method
  - [x] QR Algorithm
- [x] Solving Systems of Equations (SOE)
  - [x] Gaussian Elimination
  - [x] Steepest Descent
  - [x] Conjugate Gradient
- [x] Multi-Variable Iteration (MVI)
  - [x] Jacobi
  - [x] Gauss-Seidel
  - [x] Successive Relaxation
- [ ] Initial Value Problems (IVP)
  - [x] Forward Euler
  - [x] Backward Euler
  - [x] Improved/Modified Euler
  - [ ] Runge-Kutta
    - [ ] 2ⁿᵈ-Order
    - [x] 4ᵗʰ-Order
- [x] Boundary Value Problems (BVP)
  - [x] Linear Shooting Method
  - [x] Finite Difference Method
- [ ] Other Ordinary/Partial Differential Equations (O/PDE)
- [ ] Decision on appropriate naming convention of types/functions/methods
- [ ] Attributes and plot recipes for plotting approximations and respective errors by iteration

## Citing
See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).

## References
1. R. Burden L., D. Faires J., and A. Burden M., Numerical Analysis, 10th ed. 20 Channel Center Street, Boston, MA 02210, USA: Cengage, 2019.
