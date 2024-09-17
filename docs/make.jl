using Pkg; Pkg.precompile()
using Documenter
using NumericalMethods

DocMeta.setdocmeta!(NumericalMethods, :DocTestSetup, :(using NumericalMethods); recursive=true)

makedocs(;
    modules=[NumericalMethods],
    authors="Joby M. Anthony III",
    repo="https://github.com/jmanthony3/NumericalMethods.jl/blob/{commit}{path}#{line}",
    sitename="NumericalMethods.jl",
    doctest=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jmanthony3.github.io/NumericalMethods.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Base" => "base.md",
        "Exception" => "exception.md",
        "Single-Variable Iteration" => "singlevariableiteration.md",
        "Interpolation" => "interpolation.md",
        "Numerical Differentiation" => "numericaldifferentiation.md",
        "Numerical Integration" => "numericalintegration.md",
        "Initial-Value Problem" => "initialvalueproblem.md",
        "Multi-Variable Iteration" => "multivariableiteration.md",
        "System of Equation" => "systemofequation.md",
        "Eigenvalue" => "eigenvalue.md",
        "Boundary-Value Problem" => "boundaryvalueproblem.md"
    ],
)

deploydocs(;
    repo="github.com/jmanthony3/NumericalMethods.jl",
    devbranch="main",
)
