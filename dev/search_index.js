var documenterSearchIndex = {"docs":
[{"location":"private/#Private-Documentation","page":"Private Documentation","title":"Private Documentation","text":"","category":"section"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"Documentation for VortexFilaments.jl's private interface.","category":"page"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"See the Internals section of the manual for internal package docs covering all submodules.","category":"page"},{"location":"private/#Contents","page":"Private Documentation","title":"Contents","text":"","category":"section"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"Pages = [\"private.md\"]","category":"page"},{"location":"private/#Index","page":"Private Documentation","title":"Index","text":"","category":"section"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"Pages = [\"private.md\"]","category":"page"},{"location":"private/#Private-Interface","page":"Private Documentation","title":"Private Interface","text":"","category":"section"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"Modules = [VortexFilaments]\nPublic = false","category":"page"},{"location":"private/#Base.isinf-Tuple{StaticArrays.SArray{Tuple{2},AbstractArray,1,2}}","page":"Private Documentation","title":"Base.isinf","text":"isinf(s::StaticArrays.SArray{Tuple{2},AbstractArray,1,2}) -> Any\n\n\nChecks if s is an infinite segment.\n\n\n\n\n\n","category":"method"},{"location":"private/#Base.isinf-Tuple{VortexFilament}","page":"Private Documentation","title":"Base.isinf","text":"isinf(vf::VortexFilament) -> Bool\n\n\nChecks if vf is an infinite vortex filament.\n\n\n\n\n\n","category":"method"},{"location":"private/#VortexFilaments.infdir-Tuple{StaticArrays.SArray{Tuple{2},AbstractArray,1,2}}","page":"Private Documentation","title":"VortexFilaments.infdir","text":"infdir(s::StaticArrays.SArray{Tuple{2},AbstractArray,1,2}) -> Any\n\n\nReturns the index of the coordinate for which the first and/or second vertex of the segment s has an infinite value if it has one.\n\n\n\n\n\n","category":"method"},{"location":"private/#VortexFilaments.infdir-Tuple{VortexFilament}","page":"Private Documentation","title":"VortexFilaments.infdir","text":"infdir(vf::VortexFilament) -> Any\n\n\nReturns the index of the coordinate for which the first and/or last vertex of the vortex filament vf has an infinite value if it has one.\n\n\n\n\n\n","category":"method"},{"location":"public/#Public-Documentation","page":"Public Documentation","title":"Public Documentation","text":"","category":"section"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"Documentation for VortexFilaments.jl's public interface.","category":"page"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"See the Internals section of the manual for internal package docs covering all submodules.","category":"page"},{"location":"public/#Contents","page":"Public Documentation","title":"Contents","text":"","category":"section"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"Pages = [\"public.md\"]","category":"page"},{"location":"public/#Index","page":"Public Documentation","title":"Index","text":"","category":"section"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"Pages = [\"public.md\"]","category":"page"},{"location":"public/#Public-Interface","page":"Public Documentation","title":"Public Interface","text":"","category":"section"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"Modules = [VortexFilaments]\nPrivate = false","category":"page"},{"location":"public/#VortexFilaments.VortexFilament","page":"Public Documentation","title":"VortexFilaments.VortexFilament","text":"struct VortexFilament{T<:(AbstractArray{T,1} where T)}\n\nDefines a vortex filament of strength Γ, discretized by vertices in vertices and segments in segments, which connect the vertices. Each segment in segments should connect two subsequent entries of vertices. A vertex can be marked as a bound (e.g. to a wing) if its index appears in boundidx. Similarly, a vortex is marked as a free vertex if its index appears in freeidx.\n\nFields\n\nΓ::Real\nΓ: Strength of the vortex filament.\nvertices::Array{T,1} where T<:(AbstractArray{T,1} where T)\nvertices: Array of vertices of the vortex filament.\nsegments::Array{StaticArrays.SArray{Tuple{2},AbstractArray,1,2},1}\nsegments: Array of segments of the vortex filament.\nfreeidx::Array{Int64,1}\nfreeidx: Array of indices of the free vertices of the vortex filament.\nboundidx::Array{Int64,1}\nboundidx: Array of indices of the bound vertices of the vortex filament.\n\n\n\n\n\n","category":"type"},{"location":"public/#VortexFilaments.inducevelocity-Tuple{StaticArrays.SArray{Tuple{2},AbstractArray,1,2},Any}","page":"Public Documentation","title":"VortexFilaments.inducevelocity","text":"inducevelocity(s::StaticArrays.SArray{Tuple{2},AbstractArray,1,2}, xeval::Any) -> Any\n\n\nComputes the induced velocity by one segment s of a unit-strength vortex filament on the evaluation point xeval. In case both vertices of the segment have an infinite coordinate (which have to be in the same direction), then the non-infinite coordinates are averaged over the two vertices to calculate the distance from xeval to the segment.\n\n\n\n\n\n","category":"method"},{"location":"public/#VortexFilaments.inducevelocity-Tuple{VortexFilament,Any}","page":"Public Documentation","title":"VortexFilaments.inducevelocity","text":"inducevelocity(vf::VortexFilament, xeval::Any) -> Any\n\n\nComputes the induced velocity by the vortex filament vf on the evaluation point xeval.\n\n\n\n\n\n","category":"method"},{"location":"public/#VortexFilaments.issemiinf-Tuple{StaticArrays.SArray{Tuple{2},AbstractArray,1,2}}","page":"Public Documentation","title":"VortexFilaments.issemiinf","text":"issemiinf(s::StaticArrays.SArray{Tuple{2},AbstractArray,1,2}) -> Any\n\n\nChecks if s is a semi-infinite segment.\n\n\n\n\n\n","category":"method"},{"location":"public/#VortexFilaments.issemiinf-Tuple{VortexFilament}","page":"Public Documentation","title":"VortexFilaments.issemiinf","text":"issemiinf(vf::VortexFilament) -> Bool\n\n\nChecks if vf is a semi-infinite vortex filament.\n\n\n\n\n\n","category":"method"},{"location":"#VortexFilaments.jl","page":"Home","title":"VortexFilaments.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides tools to create and plot vortex filaments and to compute the velocity they induce in three dimensions with support for infinite and semi-infinite vortex filaments.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using VortexFilaments\nusing Plots","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package introduces the VortexFilament type, which represents a vortex filament that is discretized with vertices and segments connecting those vertices. A vortex filament can be created by calling the provided constructor,","category":"page"},{"location":"","page":"Home","title":"Home","text":"vertices = [[0.0,0.0,0.0], [0.0,1.0,0.0], [1.0,1.0,0.0], [1.0,0.0,0.0]]\nΓ = 1.0 # strength of the vortex filament\nvf = VortexFilament(Γ,vertices)","category":"page"},{"location":"","page":"Home","title":"Home","text":"which can then be plotted with the provided type recipe.","category":"page"},{"location":"","page":"Home","title":"Home","text":"plot(vf)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The velocity that the vortex filament vf induces at a location x can be computed using as inducevelocity(vf,x), which returns a 3-element vector representing the velocity vector.","category":"page"},{"location":"","page":"Home","title":"Home","text":"x = [0.5,0.5,0.5]\ninducevelocity(vf,x)","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you don't want the filament to be closed, provide the constructor with the keyword isclosed=false.","category":"page"},{"location":"","page":"Home","title":"Home","text":"vf = VortexFilament(Γ,vertices,isclosed=false)\nplot(vf)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The vortex filament can also be an infinite vortex filament or a semi-infinite vortex element. If you want to plot these filaments, you have to provide the plot axis limits for the direction in which the vortex filament extends to infinity.","category":"page"},{"location":"","page":"Home","title":"Home","text":"vertices = [[-Inf,0.0,0.0], [Inf,0.0,0.0]]\nΓ = 1.0 # strength of the vortex filament\nvf = VortexFilament(Γ,vertices) # infinite vortex filament","category":"page"},{"location":"","page":"Home","title":"Home","text":"plot(vf,xlims=[-2,2],ylims=[-2,2])","category":"page"},{"location":"","page":"Home","title":"Home","text":"vertices = [[0,0.0,0.0], [Inf,0.0,0.0]]\nΓ = 1.0 # strength of the vortex filament\nvf = VortexFilament(Γ,vertices) # semi-infinite vortex filament","category":"page"},{"location":"","page":"Home","title":"Home","text":"plot(vf,xlims=[-2,2],ylims=[-2,2])","category":"page"},{"location":"","page":"Home","title":"Home","text":"These filaments also work with the inducevelocity method. This provides the possibility to model a horseshoe vortex.","category":"page"},{"location":"","page":"Home","title":"Home","text":"b = 1\nΓ = -1.0 # sign depends on the order of the vertices\nv1 = [Inf,-b/2,0]\nv2 = [0,-b/2,0]\nv3 = [0,b/2,0]\nv4 = [Inf,b/2,0]\nvf = VortexFilament(Γ,[v1,v2,v3,v4])\nyrange1 = range(-b/2,b/2,length=20)\nxevals = [[0.0,y,0.0] for y in yrange1[2:end-1]];\nw = inducevelocity.(Ref(vf),xevals);\nnothing #hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"We will compare the induced velocity with the formula for the downwash for a horseshoe vortex.","category":"page"},{"location":"","page":"Home","title":"Home","text":"yrange2 = range(-b/2,b/2,length=100)\ndownwash(Γ,b,y) = -Γ/(4π)*b/((b/2)^2-y^2);\nnothing #hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"wvec = [[xevals[i],xevals[i]+w[i]] for i in 1:length(w)];\np = plot(vf,xlims=[-1,4],ylims=[-0.6*b,0.6*b],label=false)\nfor i in 1:length(wvec)\n    plot3d!((v->v[1]).(wvec[i]),(v->v[2]).(wvec[i]),(v->v[3]).(wvec[i]),color=:black,label=false)\nend\nplot3d!(zeros(length(yrange2)),yrange2,downwash.(1.0,b,yrange2),label=\"downwash formula\")\nplot!([],[],c=:black,label=\"inducevelocity\")","category":"page"}]
}