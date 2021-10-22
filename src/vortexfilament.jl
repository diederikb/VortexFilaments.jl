using RecipesBase

import LinearAlgebra: norm, dot, cross
import StaticArrays: SVector

export VortexFilament, Segment, getfreevertices, getboundvertices, inducevelocity, issemiinf

const Segment = SVector{2,AbstractArray}
const eps = 1e-10

"""
$(TYPEDEF)

Defines a vortex filament of strength `Γ`, discretized by vertices in `vertices` and segments in `segments`, which connect the vertices. Each segment in `segments` should connect two subsequent entries of `vertices`. A vertex can be marked as a bound (e.g. to a wing) if its index appears in `boundidx`. Similarly, a vortex is marked as a free vertex if its index appears in `freeidx`.

# Fields

$(FIELDS)
"""
struct VortexFilament{TV<:AbstractVector{<:AbstractVector},TI<:AbstractVector{Int}}
    """Γ: Strength of the vortex filament."""
    Γ::Real
    """vertices: Array of vertices of the vortex filament."""
    vertices::TV
    """segments: Array of segments of the vortex filament."""
    segments::Vector{Segment}
    """freeidx: Array of indices of the free vertices of the vortex filament."""
    freeidx::TI
    """boundidx: Array of indices of the bound vertices of the vortex filament."""
    boundidx::TI
end

"""
$(SIGNATURES)

Constructs a `VortexFilament` with strength `Γ` and vertices `v`. `boundidx` can be used to provide the indices of vertices that are considered as 'bound'. The constructor then set `freeidx` of the `VortexFilament` as the complementary list of indices. If the keyword `isclosed` is `true`, the
`segments` of the `VortexFilament` will contain a segment that connects the first and last vertex in `v`, provided that these do not contain an infinite coordinate.
"""
function VortexFilament(Γ::Real, v::TV, boundidx::TI=Int64[]; isclosed::Bool=true) where {TV<:AbstractVector{<:AbstractVector}, TI<:AbstractVector{Int}}
    infbools = map(v->in(Inf,abs.(v)),v) # boolean array with the i-th element true if the i-th element of v lies at infinity

    segments = Segment.(v[1:end-1],v[2:end])

    if !any(infbools[[1,end]]) && isclosed # if there is no vertex at infinity, close the loop by adding a segment that connects the last and first vertex
        push!(segments,Segment(v[end],v[1]))
    end

    boundidx = unique(boundidx) # remove any duplicates from boundidx
    freeidx = setdiff(1:length(v),boundidx)

    return VortexFilament{TV,TI}(Γ,v,segments,freeidx,boundidx)
end

"""
$(TYPEDSIGNATURES)

Returns the free vertices of the vortex filament `vf`.
"""
function getfreevertices(vf::VortexFilament)
    v = vf.vertices[vf.freeidx]
    return v
end

"""
$(TYPEDSIGNATURES)

Returns the bound vertices of the vortex filament `vf`.
"""
function getboundvertices(vf::VortexFilament)
    v = vf.vertices[vf.boundidx]
    return v
end

"""
$(TYPEDSIGNATURES)

Checks if `vf` is an infinite vortex filament.
"""
@inline function Base.:isinf(vf::VortexFilament)
    any(isinf,vf.segments) && return true
    semiinfsegs = findall(issemiinf,vf.segments)
    return length(semiinfsegs) > 1
end

"""
$(TYPEDSIGNATURES)

Checks if `s` is an infinite segment.
"""
@inline function Base.:isinf(s::Segment)
    return any(isinf,s[1]) && any(isinf,s[2])
end

"""
$(TYPEDSIGNATURES)

Checks if `vf` is a semi-infinite vortex filament.
"""
@inline function issemiinf(vf::VortexFilament)
    semiinfsegs = findall(issemiinf,vf.segments)
    return length(semiinfsegs) == 1
end

"""
$(TYPEDSIGNATURES)

Checks if `s` is a semi-infinite segment.
"""
@inline function issemiinf(s::Segment)
    return xor(any(isinf,s[1]),any(isinf,s[2]))
end

"""
$(TYPEDSIGNATURES)

Returns the index of the coordinate for which the first and/or last vertex of the vortex filament `vf` has an infinite value if it has one.
"""
function infdir(vf::VortexFilament)
    if isinf(vf) || issemiinf(vf)
        segdirs = infdir.(vf.segments)
        uniquesegdirs = unique(filter!(dir->dir≠nothing,segdirs))
        length(uniquesegdirs) > 1 && throw(ErrorException("Inf values in different coordinates"))
        return uniquesegdirs[1]
    else
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

Returns the index of the coordinate for which the first and/or second vertex of the segment `s` has an infinite value if it has one.
"""
function infdir(s::Segment)
    !all(length.(findall.(isinf,s)) .< 2) && throw(ErrorException("Inf values in different coordinates"))
    if issemiinf(s)
        Inf in abs.(s[1]) ? infidx = 1 : infidx = 2
        return findnext(isinf,s[infidx],1)
    elseif isinf(s)
        dir1 = findnext(isinf,s[1],1)
        dir2 = findnext(isinf,s[2],1)
        dir1 == dir2 ? dir = dir1 : throw(ErrorException("Inf values in different coordinates"))
        sign(s[1][dir]) == sign(s[2][dir]) && throw(ErrorException("Infinite segment has non-opposite infinite vertices"))
        return dir
    else
        return nothing
    end
end

"""
$(TYPEDSIGNATURES)

Computes the induced velocity by one segment `s` of a unit-strength vortex filament on the evaluation point `xeval`. In case both vertices of the segment have an infinite coordinate (which have to be in the same direction), then the non-infinite coordinates are averaged over the two vertices to calculate the distance from xeval to the segment.
"""
function inducevelocity(s::Segment,xeval)
    if s[1] == s[2]
        return zeros(3)
    elseif issemiinf(s)
        #=
        In this case we calculate the induced velocity by splitting the semi-infinite segment into two segments using the point `v`, which is the closest point to xeval on the axis of the semi-infinite segment. At the end we sum the velocities induced by the two segments using two multipliers `K1` and `K2` to determine the sign of the velocities.

        The first segment `vfiniteseg` is the segment from the finite end of `s` to `v`. Because the sign of the velocity induced by this segment flips if `v` lies outside of `s` (because then this velocity has to be subtracted) and if the order of the two vertices of `s` is flipped, we multiply it by `K1` and `K2` when computing `vnet`.

        The second segment `vseminfseg` is the segment from the infinite end of `s` to `v`. Because the sign of the velocity induced by this segment flips if the order of the two vertices of `s` is flipped, we multiply it by `K1` when computing `vnet`.
        =#

        infidx = findnext(v->in(Inf,abs.(v)),s,1)
        infidx == 1 ? (finidx = 2; K1 = -1) : (finidx = 1; K1 = 1)

        dir = infdir(s)
        edir = zeros(3)
        edir[dir] = sign(s[infidx][dir])

        v,t = _closestpointonfinitesegmentaxis(Segment(s[finidx],s[finidx]+edir),xeval)
        t < 0 ? K2 = -1 : K2 = 1 # if t < 0, v does not lie on the segment
        vfiniteseg = inducevelocity(Segment(v,s[finidx]),xeval)

        w = xeval - v

        if norm(w) < eps
            vnet = zeros(3)
        else
            vseminfseg = -1/(4π)*cross(w,edir)/dot(w,w)
            vnet = K1*K2*vfiniteseg + K1*vseminfseg
        end
        return vnet

    elseif isinf(s)
        dir = infdir(s)
        edir = zeros(3)
        edir[dir] = sign(s[2][dir])

        centerpoint = (s[1]+s[2])/2
        centerpoint[dir] = 0.0

        v,t = _closestpointonfinitesegmentaxis(Segment(centerpoint,centerpoint+edir),xeval)

        w = xeval - v

        if norm(w) < eps
            vinfseg = zeros(3)
        else
            vinfseg = -1/(2π)*cross(w,edir)/dot(w,w)
        end
        return vinfseg
    else
        r0 = s[2] - s[1]
        r1 = xeval - s[1]
        r2 = xeval - s[2]
        r1norm = norm(r1)
        r2norm = norm(r2)

        r1r2 = cross(r1,r2)
        r1r2sq = dot(r1r2,r1r2)

        if r1norm < eps || r2norm < eps || r1r2sq < eps
            vseg = zeros(3)
        else
            r0r1 = dot(r0,r1)
            r0r2 = dot(r0,r2)
            K = 1.0/(4.0*π*r1r2sq)*(r0r1/r1norm - r0r2/r2norm)
            vseg = K*r1r2
        end

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

@inline function _closestpointonfinitesegmentaxis(s::Segment,p)
    a = s[1]
    b = s[2]
    ab = b-a
    ap = p-a
    t = dot(ap,ab)/dot(ab,ab)
    return a + t*ab, t
end

lims(dir) = (:xlims, :ylims, :zlims)[dir]

@recipe function f(vf::VortexFilament)
    segments = vf.segments
    if isinf(vf) || issemiinf(vf)
        segments = deepcopy(vf.segments)
        for s in segments
            if isinf(s) || issemiinf(s)
                inftydir = infdir(s)
                noninftydir = filter!(x->x≠inftydir,[1,2,3])
                inflims = get(plotattributes,lims(inftydir),:none)
                inflims == :none && throw(ArgumentError("The filament has a vertex at infinity, please provide axis limits"))
                if isinf(s)
                    centerpoint = 0.5*(s[1] + s[2])
                    s[1][noninftydir] .= centerpoint[noninftydir]
                    s[2][noninftydir] .= centerpoint[noninftydir]
                end
                if isinf(s[1][inftydir])
                    sign(s[1][inftydir]) == -1 ? limitidx = 1 : limitidx = 2
                    s[1][inftydir] = inflims[limitidx]
                    s[1][noninftydir] .= s[2][noninftydir]
                end
                if isinf(s[2][inftydir])
                    sign(s[2][inftydir]) == -1 ? limitidx = 1 : limitidx = 2
                    s[2][inftydir] = inflims[limitidx]
                    s[2][noninftydir] .= s[1][noninftydir]
                end
            end
        end
    end
    vertices = vcat((s->s[1]).(segments),[segments[end][end]])
    return (v->v[1]).(vertices), (v->v[2]).(vertices), (v->v[3]).(vertices)
end

Base.:(==)(a::VortexFilament, b::VortexFilament) = isequal(a.Γ, b.Γ) && isequal(a.vertices, b.vertices) &&  isequal(a.segments, b.segments) && isequal(a.freeidx, b.freeidx) && isequal(a.boundidx, b.boundidx)

function Base.show(io::IO, vf::VortexFilament)
    isinf(vf) ? type = "Infinite" : issemiinf(vf) ? type = "Semi-infinite" : type = "Finite"
    nv = length(vf.vertices)
    ns = length(vf.segments)

    println(io, "$type vortex filament with $((nv == 1) ? "1 vertex" : "$nv vertices"), $((ns == 1) ? "1 segment" : "$ns segments"), and strength Γ = $(vf.Γ)")
end
