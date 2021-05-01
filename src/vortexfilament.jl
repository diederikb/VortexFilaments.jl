using RecipesBase

import StaticArrays: SVector
import LinearAlgebra: norm, dot, cross

export VortexFilament, Segment, getfreevertices, inducevelocity

const sigma = 0.005;
const sigsq = sigma*sigma;
const Segment = SVector{2,Vertex}

"""
$(TYPEDEF)

Defines a vortex filament of strength `Γ` consisting segments that connect the N vertices in `vertices`. Each segment in `segments` connects two subsequent entries of `vertices`. In case the first or last vertex has a coordinate value at infinity (`Inf`), there will be N-1 segments. Otherwise, there will be an N-th segment connecting the last and first vertex and thereby closing the filament. A vertex can be marked as a bound to a solid (e.g. a wing) if its index appears in `boundidx`. Similarly, a vortex is marked as a free vertex if its index appears in `freeidx`. If either the first or last vertex lies at infinity and the other does not, that other vertex is a bound vertex.

# Fields

$(TYPEDFIELDS)
"""
struct VortexFilament
    """Γ: Strength of the vortex filament."""
    Γ::Real
    """vertices: Array of vertices of the vortex filament."""
    vertices::Vector{AbstractVector}
    """segments: Array of segments of the vortex filament."""
    segments::Vector{Segment}
    """freeidx: Array of indices of the free vertices of the vortex filament."""
    freeidx::Vector{Int}
    """boundidx: Array of indices of the bound vertices of the vortex filament."""
    boundidx::Vector{Int}
end

    function VortexFilament(Γ::Real, vertices::Vector{<:AbstractVector}, boundidx::Vector{Int}=Int64[])
    infbools = map(v->in(Inf,abs.(v)),vertices) # boolean array with the i-th element true if the i-th element of vertices lies at infinity
    @assert length(vertices) > 1 "The vortex filament has to contain at least two vertices."
    @assert !in(true,infbools[2:end-1]) "Only the first or last vertex can have an infinite coordinate value."
    @assert !(in(1,boundidx) && infbools[1]) "The first vertex lies at infinity and can therefore not be a bound vertex."
    @assert !(in(length(vertices),boundidx) && infbools[end]) "The last vertex lies at infinity and can therefore not be a bound vertex."

    if vertices[1] == vertices[end] # if the first and last vertex in vertices are equal, remove the last one
        pop!(vertices)
    end

    segments = Segment.(vertices[1:end-1],vertices[2:end])

    if !in(true,infbools) # if there is no vertex at infinity, close the loop by adding a segment that connects the last and first vertex
        push!(segments,Segment(vertices[end],vertices[1]))
    end

    if infbools[1] && !infbools[end] # if either the first or last vertex lies at infinity and the other does not, the other vertex has to be a bound vertex
        push!(boundidx,length(vertices))
    elseif !infbools[1] && infbools[end]
        push!(boundidx,1)
    end

    boundidx = unique(boundidx) # remove any duplicates from boundidx
    freeidx = setdiff(1:length(vertices),boundidx)

    return VortexFilament(Γ,vertices,segments,freeidx,boundidx)
end

"""
$(TYPEDSIGNATURES)

Returns the vertices of the vortex filament `vf` that are free.
"""
function getfreevertices(vf::VortexFilament)
    vf.vertices[vf.freeidx]
end

"""
$(TYPEDSIGNATURES)

Computes the induced velocity by one segment `s` of a unit-strength vortex filament on the evaluation point `xeval`.
"""
function inducevelocity(s::Segment,xeval)
    r = Ref(xeval) .- s;
    xnorm = sqrt.(dot.(r,r) .+ sigsq);
    xdir = r./xnorm;

    x21  = s[2] - s[1];
    lensq = dot(x21,x21);

    # b = rj x rj+1
    b = cross(r[1],r[2]);
    bsq = dot(b,b)+sigsq*lensq;

    # f = rj/|rj| - rj+1/|rj+1|
    f = xdir[1] - xdir[2];

    g = -dot(x21,f);
    boverbsq = b/bsq;

    # Velocity
    vseg = -boverbsq*g;
    return vseg
end

"""
$(TYPEDSIGNATURES)

Computes the induced velocity by the vortex filament `vf` on the evaluation point `xeval`.
"""
function inducevelocity(vf::VortexFilament,xeval)
    vel = sum(inducevelocity.(vf.segments,Ref(xeval)))
    vel *= vf.Γ/(4*pi)
end

@recipe f(vf::VortexFilament) = map(p->p.x, vcat(vf.vertices,[vf.vertices[1]])), map(p->p.y, vcat(vf.vertices,[vf.vertices[1]])), map(p->p.z, vcat(vf.vertices,[vf.vertices[1]]))

Base.:(==)(a::VortexFilament, b::VortexFilament) = isequal(a.Γ, b.Γ) && isequal(a.vertices, b.vertices) &&  isequal(a.segments, b.segments) && isequal(a.freeidx, b.freeidx) && isequal(a.boundidx, b.boundidx)
