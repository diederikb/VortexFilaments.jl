import StaticArrays: FieldVector

export Vertex, add!

"""
$(TYPEDEF)

Defines a vertex with at position (`x`,`y`,`z`).

# Fields

$(TYPEDFIELDS)
"""
mutable struct Vertex <: FieldVector{3, Float64}
    x::Float64
    y::Float64
    z::Float64
end

"""
$(TYPEDSIGNATURES)

Adds a vector `Δv` to the position of vertex `v`.
"""
function add!(v::Vertex,Δv)
    v .+= Δv
end
