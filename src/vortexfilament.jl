using RecipesBase

import StaticArrays: SVector
import LinearAlgebra: norm, dot, cross

export VortexFilament, Segment, getfreevertices, inducevelocity

const sigma = 0.005;
const sigsq = sigma*sigma;
const Segment = SVector{2,Vertex}

"""
$(TYPEDEF)

Defines a vortex filament of strength `Γ` consisting of N segments that connect the N vertices in `vertices`. Each segment in `segments` connects two subsequent entries of `vertices`, with the last segment connecting the last and first entry and thereby closing the filament. Vertices with indices in `freeidx` are free and can move with the flow. Vertices with indices in `boundidx` are bound to a wing.

# Fields

$(TYPEDFIELDS)
"""
struct VortexFilament
    """Γ: Strength of the vortex filament."""
    Γ::Real
    """vertices: Array of vertices of the vortex filament."""
    vertices::Vector{Vertex}
    """segments: Array of segments of the vortex filament."""
    segments::Vector{Segment}
    """freeidx: Array of indices of the free vertices of the vortex filament."""
    freeidx::Vector{Int}
    """boundidx: Array of indices of the bound vertices of the vortex filament."""
    boundidx::Vector{Int}
end

function VortexFilament(Γ::Real, vertices::Vector{Vertex}, boundidx::Vector=[])
    freeidx = setdiff(1:length(vertices),boundidx)
    segments = Segment.(vertices,vcat(vertices[2:end],[vertices[1]]))
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

Returns the induced velocity by one segment `s` of a unit-strength vortex filament on the evaluation point `xeval`.
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

Returns the induced velocity by the vortex filament `vf` on the evaluation point `xeval`.
"""
function inducevelocity(vf::VortexFilament,xeval)
    vel = sum(inducevelocity.(vf.segments,Ref(xeval)))
    vel *= vf.Γ/(4*pi)
end

@recipe f(vf::VortexFilament) = map(p->p.x, vcat(vf.vertices,[vf.vertices[1]])), map(p->p.y, vcat(vf.vertices,[vf.vertices[1]])), map(p->p.z, vcat(vf.vertices,[vf.vertices[1]]))
