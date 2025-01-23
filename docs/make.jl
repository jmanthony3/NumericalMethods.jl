using Pkg; Pkg.precompile()
using Documenter
using LUSE_ENGR701_704_NumericalMethods

DocMeta.setdocmeta!(LUSE_ENGR701_704_NumericalMethods, :DocTestSetup, :(using LUSE_ENGR701_704_NumericalMethods); recursive=true)

makedocs(;
    modules=[LUSE_ENGR701_704_NumericalMethods],
    authors="Joby M. Anthony III",
    repo="https://github.com/jmanthony3/LUSE_ENGR701_704_NumericalMethods.jl/blob/{commit}{path}#{line}",
    sitename="LUSE_ENGR701_704_NumericalMethods.jl",
    doctest=false,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jmanthony3.github.io/LUSE_ENGR701_704_NumericalMethods.jl",
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
    repo="github.com/jmanthony3/LUSE_ENGR701_704_NumericalMethods.jl",
    devbranch="main",
)
