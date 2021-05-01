using VortexFilaments
using Test

const GROUP = get(ENV, "GROUP", "All")

notebookdir = "../examples"

@testset "Vortex filament" begin
    include("vortexfilament.jl")
end
