using Documenter
using VortexFilaments

ENV["GKSwstype"] = "nul"

makedocs(
    modules = [VortexFilaments],
    sitename = "VortexFilaments.jl",
    doctest = true,
    clean = true,
    pages = [
        "Home" => "index.md",
        "Library" => ["public.md",
                      "private.md"]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax2(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    )
)

deploydocs(
    repo = "github.com/diederikb/VortexFilaments.jl.git",
    target = "build"
)
