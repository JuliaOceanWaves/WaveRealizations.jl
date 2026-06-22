using Documenter
using WaveRealizations

makedocs(
    sitename = "WaveRealizations.jl",
    modules = [WaveRealizations],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaOceanWaves/WaveRealizations.jl.git",
    push_preview = true,
    versions = ["latest" => "v^", "v#"])
