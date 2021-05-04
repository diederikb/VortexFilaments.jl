using RecipesBase

import StaticArrays: SVector
import LinearAlgebra: norm, dot, cross

export VortexFilament, Segment, getfreevertices, inducevelocity

const sigma = 0.005;
const sigsq = sigma*sigma;
const Segment = SVector{2,AbstractArray}

"""
$(TYPEDEF)

Defines a vortex filament of strength `Γ` consisting segments that connect the N vertices in `vertices`. Each segment in `segments` connects two subsequent entries of `vertices`. In case the first or last vertex has a coordinate value at infinity (`Inf`), there will be N-1 segments. Otherwise, there will be an N-th segment connecting the last and first vertex and thereby closing the filament. A vertex can be marked as a bound (e.g. to a wing) if its index appears in `boundidx`. Similarly, a vortex is marked as a free vertex if its index appears in `freeidx`. If either the first or last vertex lies at infinity and the other does not, that other vertex is a bound vertex.

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
    @assert all(length.(vertices) .== 3) "Each vertex should have three coordinates."
    @assert all(length.(findall.(isinf,vertices)) .< 2) "A vertex cannot have more than one infinite coordinate value."
    @assert !in(true,infbools[2:end-1]) "Only the first or last vertex can have an infinite coordinate value."
    @assert !(in(1,boundidx) && infbools[1]) "The first vertex lies at infinity and can therefore not be a bound vertex."
    @assert !(in(length(vertices),boundidx) && infbools[end]) "The last vertex lies at infinity and can therefore not be a bound vertex."
    if length(vertices) == 2 && all(infbools)
        dir1 = findnext(isinf,vertices[1],1)
        dir2 = findnext(isinf,vertices[2],1)
        @assert (dir1 == dir2 && vertices[1][dir1] == -vertices[2][dir2]) "A filament with just two vertices, both at infinity, should have its vertices at minus and plus infinity in the same direction."
    end

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
    if s[1] == s[2]
        return zeros(3)
    elseif xor(in(Inf,abs.(s[1])),in(Inf,abs.(s[2])))
        #=
        In this case we calculate the induced velocity by splitting the semi-infinite segment into two segments using the point `v`, which is the closest point to xeval on the axis of the semi-infinite segment. At the end we sum the velocities induced by the two segments using two multipliers `K1` and `K2` to determine the sign of the velocities.

        The first segment `vfiniteseg` is the segment from the finite end of `s` to `v`. Because the sign of the velocity induced by this segment flips if `v` lies outside of `s` (because then this velocity has to be subtracted) and if the order of the two vertices of `s` is flipped, we multiply it by `K1` and `K2` when computing `vnet`.

        The second segment `vseminfseg` is the segment from the infinite end of `s` to `v`. Because the sign of the velocity induced by this segment flips if the order of the two vertices of `s` is flipped, we multiply it by `K1` when computing `vnet`.
        =#

        infidx = findnext(v->in(Inf,abs.(v)),s,1)
        infidx == 1 ? (finidx = 2; K1 = -1) : (finidx = 1; K1 = 1)

        dir = findnext(isinf,s[infidx],1)
        edir = zeros(3)
        edir[dir] = sign(s[infidx][dir])

        v,t = _closestpointonfinitesegmentaxis(Segment(s[finidx],s[finidx]+edir),xeval)
        t < 0 ? K2 = -1 : K2 = 1 # if t < 0, v does not lie on the segment
        vfiniteseg = inducevelocity(Segment(v,s[finidx]),xeval)

        vseminfseg = 1/(4π)*cross(xeval-v,edir)/(dot(xeval-v,xeval-v)+ sigsq)
        vnet = K1*K2*vfiniteseg + K1*vseminfseg
        return vnet

    elseif in(Inf,abs.(s[1])) && in(Inf,abs.(s[2]))
        dir = findnext(isinf,s[1],1)
        edir = zeros(3)
        edir[dir] = sign(s[2][dir])

        centerpoint = (s[1]+s[2])/2
        centerpoint[dir] = 0.0

        v,t = _closestpointonfinitesegmentaxis(Segment(centerpoint,centerpoint+edir),xeval)
        vinfseg = 1/(2π)*cross(xeval-v,edir)/(dot(xeval-v,xeval-v) + sigsq)
        return vinfseg
    else
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
        vseg = -1/(4π)*boverbsq*g;

        return vseg
    end
end

"""
$(TYPEDSIGNATURES)

Computes the induced velocity by the vortex filament `vf` on the evaluation point `xeval`.
"""
function inducevelocity(vf::VortexFilament,xeval)
    vel = sum(inducevelocity.(vf.segments,Ref(xeval)))
    vel *= vf.Γ
end

function _closestpointonsemiinfinitesegment(s::Segment,p)
    infidx = findnext(v->in(Inf,abs.(v)),s,1)
    infidx == 1 ? finidx = 2 : finidx = 1
    dir = findnext(isinf,s[infidx],1)
    a = s[finidx]
    segdir = zeros(3)
    segdir[dir] = sign(s[infidx][dir])
    ap = p-a
    t = dot(ap,segdir)
    if t < 0
        return a # p lies further than the segment, so the finite endpoint of the segment is the closest point.
    else
        return a + t*segdir
    end
end

function _closestpointonfinitesegmentaxis(s::Segment,p)
    a = s[1]
    b = s[2]
    ab = b-a
    ap = p-a
    t = dot(ap,ab)/dot(ab,ab)
    return a + t*ab, t
end

@recipe f(vf::VortexFilament) = map(p->p.x, vcat(vf.vertices,[vf.vertices[1]])), map(p->p.y, vcat(vf.vertices,[vf.vertices[1]])), map(p->p.z, vcat(vf.vertices,[vf.vertices[1]]))

Base.:(==)(a::VortexFilament, b::VortexFilament) = isequal(a.Γ, b.Γ) && isequal(a.vertices, b.vertices) &&  isequal(a.segments, b.segments) && isequal(a.freeidx, b.freeidx) && isequal(a.boundidx, b.boundidx)
