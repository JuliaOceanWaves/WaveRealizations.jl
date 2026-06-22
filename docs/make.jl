using Documenter
using WaveRealizations

DocMeta.setdocmeta!(
    WaveRealizations,
    :DocTestSetup,
    :(using WaveRealizations);
    recursive = true
)

makedocs(
    sitename = "WaveRealizations.jl",
    modules = [WaveRealizations],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaOceanWaves/WaveRealizations.jl.git",
    devbranch = "main",
    push_preview = true
)
