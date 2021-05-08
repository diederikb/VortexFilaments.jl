using VortexFilaments
using Test

const GROUP = get(ENV, "GROUP", "All")

ENV["GKSwstype"] = "nul"

include("vortexfilament.jl")
