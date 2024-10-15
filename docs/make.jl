using SteadyWaves
using Documenter
using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

DocMeta.setdocmeta!(SteadyWaves, :DocTestSetup, :(using SteadyWaves); recursive=true)

makedocs(;
    modules=[SteadyWaves],
    authors="Maciej Paprota",
    sitename="SteadyWaves.jl",
    plugins=[bib],
    format=Documenter.HTML(;
        canonical="https://mcpaprota.github.io/SteadyWaves.jl",
        edit_link="main",
        assets=String["assets/citations.css"],
    ),
    pages=[
        "Home" => "index.md",
        "Guide" => "guide.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
)

deploydocs(;
    repo="github.com/mcpaprota/SteadyWaves.jl",
)
