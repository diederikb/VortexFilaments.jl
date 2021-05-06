# VortexFilaments.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://diederikb.github.io/VortexFilaments.jl/dev)
[![Build Status](https://github.com/diederikb/VortexFilaments.jl/workflows/CI/badge.svg)](https://github.com/diederikb/VortexFilaments.jl/actions)
[![Coverage](https://codecov.io/gh/diederikb/VortexFilaments.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/diederikb/VortexFilaments.jl)

This package provides tools to create and plot vortex filaments and to compute the velocity they induce in three dimensions with support for infinite and semi-infinite vortex filaments.

<!-- **VortexFilaments.jl** is registered in the general Julia registry. To install, type
e.g.,
```julia
] add VortexFilaments
```

Then, in any version, type
```julia
julia> using VortexFilaments
``` -->

The package introduces the `VortexFilament` type, which represents a vortex filament that is discretized with vertices and segments connecting those vertices. A vortex filament can be created by calling the provided constructor,

```julia
vertices = [[0.0,0.0,0.0], [0.0,1.0,0.0], [1.0,1.0,0.0], [1.0,0.0,0.0]]
Γ = 1.0 # strength of the vortex filament
vf = VortexFilament(Γ,vertices)
```

which can then be plotted with the provided type recipe.

```julia
plot(vf)
```
![closedfinitevl](https://user-images.githubusercontent.com/26737762/117224622-4dc74680-adc5-11eb-97fa-fdc9a33779c6.png)

The velocity that the vortex filament `vf` induces at a location `x` can be computed using as `inducevelocity(vf,x)`, which returns a 3-element vector representing the velocity vector.
```julia
x = [0.5,0.5,0.5]
inducevelocity(vf,x)
```

If you don't want the filament to be closed, provide the constructor with the keyword `isclosed=false`.

```julia
vf = VortexFilament(Γ,vertices,isclosed=false)
plot(vf)
```
![openfinitevl](https://user-images.githubusercontent.com/26737762/117224654-60418000-adc5-11eb-9b76-c6fe3e6de4bb.png)

The vortex filament can also be an infinite vortex filament or a semi-infinite vortex element. If you want to plot these filaments, you have to provide the plot axis limits for the direction in which the vortex filament extends to infinity.

```julia
vertices = [[-Inf,0.0,0.0], [Inf,0.0,0.0]]
Γ = 1.0 # strength of the vortex filament
vf = VortexFilament(Γ,vertices) # infinite vortex filament
plot(vf,xlims=[-2,2],ylims=[-2,2])
```
![infinitevl](https://user-images.githubusercontent.com/26737762/117224672-6fc0c900-adc5-11eb-9fbf-065aa62e1519.png)

```julia
vertices = [[0,0.0,0.0], [Inf,0.0,0.0]]
Γ = 1.0 # strength of the vortex filament
vf = VortexFilament(Γ,vertices) # semi-infinite vortex filament
plot(vf,xlims=[-2,2],ylims=[-2,2])
```
![semiinfinitevl](https://user-images.githubusercontent.com/26737762/117224686-75b6aa00-adc5-11eb-8707-3d0a258b6c3a.png)

These filaments also work with the `inducevelocity` method. This provides the possibility to model a horseshoe vortex.
```julia
b = 1
Γ = -1.0 # sign depends on the order of the vertices
v1 = [Inf,-b/2,0]
v2 = [0,-b/2,0]
v3 = [0,b/2,0]
v4 = [Inf,b/2,0]
vf = VortexFilament(Γ,[v1,v2,v3,v4])
yrange1 = range(-b/2,b/2,length=20)
xevals = [[0.0,y,0.0] for y in yrange1[2:end-1]];
w = inducevelocity.(Ref(vf),xevals);
```

We will compare the induced velocity with the formula for the downwash for a horseshoe vortex.
```julia
yrange2 = range(-b/2,b/2,length=100)
downwash(Γ,b,y) = -Γ/(4π)*b/((b/2)^2-y^2);
```

```julia
wvec = [[xevals[i],xevals[i]+w[i]] for i in 1:length(w)];
p = plot(vf,xlims=[-1,4],ylims=[-0.6*b,0.6*b],label=false)
for i in 1:length(wvec)
    plot3d!((v->v[1]).(wvec[i]),(v->v[2]).(wvec[i]),(v->v[3]).(wvec[i]),color=:black,label=false)
end
plot3d!(zeros(length(yrange2)),yrange2,downwash.(1.0,b,yrange2),label="downwash formula")
plot!([],[],c=:black,label="inducevelocity")
```
![horseshoe](https://user-images.githubusercontent.com/26737762/117225243-b95de380-adc6-11eb-85d1-5f213e04dc6f.png)
