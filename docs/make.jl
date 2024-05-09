using Documenter
using SynthDiag

makedocs(;
    modules=[SynthDiag],
    format=Documenter.HTML(),
    sitename="SynthDiag",
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/ProjectTorreyPines/SynthDiag.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="dev",
    versions=["stable" => "v^", "v#.#"],
)
