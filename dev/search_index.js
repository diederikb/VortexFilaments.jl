var documenterSearchIndex = {"docs":
[{"location":"private/#Private-Documentation","page":"Private Documentation","title":"Private Documentation","text":"","category":"section"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"Documentation for VortexFilaments.jl's private interface.","category":"page"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"See the Internals section of the manual for internal package docs covering all submodules.","category":"page"},{"location":"private/#Contents","page":"Private Documentation","title":"Contents","text":"","category":"section"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"Pages = [\"private.md\"]","category":"page"},{"location":"private/#Index","page":"Private Documentation","title":"Index","text":"","category":"section"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"Pages = [\"private.md\"]","category":"page"},{"location":"private/#Private-Interface","page":"Private Documentation","title":"Private Interface","text":"","category":"section"},{"location":"private/","page":"Private Documentation","title":"Private Documentation","text":"Modules = [VortexFilaments]\nPublic = false","category":"page"},{"location":"public/#Public-Documentation","page":"Public Documentation","title":"Public Documentation","text":"","category":"section"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"Documentation for VortexFilaments.jl's public interface.","category":"page"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"See the Internals section of the manual for internal package docs covering all submodules.","category":"page"},{"location":"public/#Contents","page":"Public Documentation","title":"Contents","text":"","category":"section"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"Pages = [\"public.md\"]","category":"page"},{"location":"public/#Index","page":"Public Documentation","title":"Index","text":"","category":"section"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"Pages = [\"public.md\"]","category":"page"},{"location":"public/#Public-Interface","page":"Public Documentation","title":"Public Interface","text":"","category":"section"},{"location":"public/","page":"Public Documentation","title":"Public Documentation","text":"Modules = [VortexFilaments]\nPrivate = false","category":"page"},{"location":"public/#VortexFilaments.VortexFilament","page":"Public Documentation","title":"VortexFilaments.VortexFilament","text":"struct VortexFilament\n\nDefines a vortex filament of strength Γ consisting segments that connect the N vertices in vertices. Each segment in segments connects two subsequent entries of vertices. In case the first or last vertex has a coordinate value at infinity (Inf), there will be N-1 segments. Otherwise, there will be an N-th segment connecting the last and first vertex and thereby closing the filament. A vertex can be marked as a bound to a solid (e.g. a wing) if its index appears in boundidx. Similarly, a vortex is marked as a free vertex if its index appears in freeidx. If either the first or last vertex lies at infinity and the other does not, that other vertex is a bound vertex.\n\nFields\n\nΓ::Real\nΓ: Strength of the vortex filament.\nvertices::Array{AbstractArray{T,1} where T,1}\nvertices: Array of vertices of the vortex filament.\nsegments::Array{StaticArrays.SArray{Tuple{2},AbstractArray,1,2},1}\nsegments: Array of segments of the vortex filament.\nfreeidx::Array{Int64,1}\nfreeidx: Array of indices of the free vertices of the vortex filament.\nboundidx::Array{Int64,1}\nboundidx: Array of indices of the bound vertices of the vortex filament.\n\n\n\n\n\n","category":"type"},{"location":"public/#VortexFilaments.getfreevertices-Tuple{VortexFilament}","page":"Public Documentation","title":"VortexFilaments.getfreevertices","text":"getfreevertices(vf::VortexFilament) -> Array{AbstractArray{T,1} where T,1}\n\n\nReturns the vertices of the vortex filament vf that are free.\n\n\n\n\n\n","category":"method"},{"location":"public/#VortexFilaments.inducevelocity-Tuple{StaticArrays.SArray{Tuple{2},AbstractArray,1,2},Any}","page":"Public Documentation","title":"VortexFilaments.inducevelocity","text":"inducevelocity(s::StaticArrays.SArray{Tuple{2},AbstractArray,1,2}, xeval::Any) -> Any\n\n\nComputes the induced velocity by one segment s of a unit-strength vortex filament on the evaluation point xeval.\n\n\n\n\n\n","category":"method"},{"location":"public/#VortexFilaments.inducevelocity-Tuple{VortexFilament,Any}","page":"Public Documentation","title":"VortexFilaments.inducevelocity","text":"inducevelocity(vf::VortexFilament, xeval::Any) -> Any\n\n\nComputes the induced velocity by the vortex filament vf on the evaluation point xeval.\n\n\n\n\n\n","category":"method"},{"location":"#LiftingLineFlow.jl","page":"Home","title":"LiftingLineFlow.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A set of tools to solve flow past finite wings using steady and unsteady lifting line theory.","category":"page"}]
}