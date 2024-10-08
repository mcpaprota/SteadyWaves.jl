using SteadyWaves
using Documenter

DocMeta.setdocmeta!(SteadyWaves, :DocTestSetup, :(using SteadyWaves); recursive=true)

makedocs(;
    modules=[SteadyWaves],
    authors="Maciej Paprota",
    sitename="SteadyWaves.jl",
    format=Documenter.HTML(;
        canonical="https://mcpaprota.github.io/SteadyWaves.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/mcpaprota/SteadyWaves.jl",
)
