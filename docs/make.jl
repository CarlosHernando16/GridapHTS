using Documenter
using GridapHTS

makedocs(
    sitename = "GridapHTS.jl",
    authors  = "Magnestar Tech",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    modules  = [GridapHTS],
    pages    = [
        "Home" => "index.md",
        "Formulations" => "formulations.md",
        "Superconductivity Modeling Context" => "superconductivity_modeling_context.md",
        "Benchmarks" => "benchmarks.md",
        "API Reference" => "api.md",
        "Tutorials" => [
            "tutorials/getting_started.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/magnestar/GridapHTS.jl.git",
    devbranch = "main",
)
