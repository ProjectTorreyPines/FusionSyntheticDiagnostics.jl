using Documenter
using FusionSyntheticDiagnostics

makedocs(;
    modules=[FusionSyntheticDiagnostics],
    format=Documenter.HTML(),
    sitename="FusionSyntheticDiagnostics",
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/ProjectTorreyPines/FusionSyntheticDiagnostics.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    versions=["stable" => "v^", "v#.#"],
)
