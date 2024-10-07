using SteadyWaves
using Documenter

DocMeta.setdocmeta!(SteadyWaves, :DocTestSetup, :(using SteadyWaves); recursive=true)

makedocs(;
    modules=[SteadyWaves],
    authors="Maciej Paprota",
    sitename="SteadyWaves.jl",
    format=Documenter.HTML(;
        canonical="https://mcpaprota.github.io/SteadyWaves.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mcpaprota/SteadyWaves.jl",
    devbranch="master",
)
