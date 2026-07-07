using Documenter
using WaveRealizations

DocMeta.setdocmeta!(SpectralSuperpositions, :DocTestSetup, :(using SpectralSuperpositions); recursive = true)
#bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

makedocs(
    sitename = "WaveRealizations.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        assets = String["src/index.css"]
    ),
    modules = [SpectralSuperpositions],
    pages = ["Home" => "index.md",
            "API" => "api.md"]    #plugins = [bib],
)

deploydocs(
    repo = "github.com/JuliaOceanWaves/WaveRealizations.jl.git",
    versions = ["stable" => "v^", "v#"],
    devbranch = "main",
    push_preview = true
)