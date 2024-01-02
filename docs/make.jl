using NumericalMethods
using Documenter

DocMeta.setdocmeta!(NumericalMethods, :DocTestSetup, :(using NumericalMethods); recursive=true)

makedocs(;
    modules=[NumericalMethods],
    authors="Joby M. Anthony III",
    repo="https://github.com/jmanthony3/NumericalMethods.jl/blob/{commit}{path}#{line}",
    sitename="NumericalMethods.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jmanthony3.github.io/NumericalMethods.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jmanthony3/NumericalMethods.jl",
    devbranch="master",
)
