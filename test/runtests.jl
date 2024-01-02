using NumericalMethods
using Test

@testset "NumericalMethods.jl" begin
    # # ENGR 705-001: Finite Element Analysis - Homework 3
    # ## 3.1
    # A, E, L = 4, 30e6, 30 # [in², psi, in]
    # elements = Dict("1"=>(1, 2), "2"=>(2, 3), "3"=>(3, 4))
    # nodeboundaryconditions = (1, 4)
    # F_applied = [(2, 5e3), (3, -10e3)]
    # preparedelements = prepareelements_line(elements, A, E, L)
    # U = solve(preparedelements, nodeboundaryconditions, F_applied).U
    # @test round.(U, digits=6) == [0, 0, -0.00125, 0]
    # ## 3.2
    # A, E, L = 4e-4, 210e9, 2 # [m², Pa, m]
    # elements = Dict("1"=>(1, 2), "2"=>(2, 3))
    # nodeboundaryconditions = [1, (3, 25e-3, 2)]
    # F_applied = [(2, -5e3)]
    # preparedelements = prepareelements_line(elements, A, E, L)
    # U = solve(preparedelements, nodeboundaryconditions, F_applied).U
    # @test round.(U, digits=6) == [0, 0.012440, 0.025]
    # ## 3.3
    # A, E, L = 1, 10e6, 100          # [in², psi, in]
    # angles = [120, 0, 210]          # [°]
    # L = L ./ abs.(cosd.(angles))
    # elements = Dict("1"=>(1, 2), "2"=>(1, 3), "3"=>(1, 4))
    # nodeboundaryconditions = [2, 3, 4]
    # F_applied = [(1, (1e3, 1e3))]   # (node, (force_x [lb], force_y [lb]))
    # preparedelements = prepareelements_line(elements, A, E, L, angles, dims=2)
    # sol = solve(preparedelements, nodeboundaryconditions, F_applied, dims=2)
    # σ = axialstresses(preparedelements, sol.U, dims=2)
    # @test round.(σ, digits=6) == [-577.350269, -422.649731, 1000]

    # Write your tests here.
end
