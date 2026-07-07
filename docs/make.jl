using Documenter
using WaveRealizations

DocMeta.setdocmeta!(WaveRealizations, :DocTestSetup, :(using WaveRealizations); recursive = true)
#bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

makedocs(
    sitename = "WaveRealizations.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        assets = String["src/index.css"]
    ),
    modules = [WaveRealizations],
    pages = ["Home" => "index.md",
            "API" => "api.md"]    #plugins = [bib],
)

deploydocs(
    repo = "github.com/JuliaOceanWaves/WaveRealizations.jl.git",
    versions = ["stable" => "v^", "v#"],
    devbranch = "main",
    push_preview = true
)