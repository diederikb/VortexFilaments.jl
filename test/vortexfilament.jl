using Random
using LinearAlgebra: norm

@testset "VortexFilament constructor" begin
    v1 = rand(3)
    v2 = rand(3)
    v3 = rand(3)

    @test VortexFilament(1.0,[v1,v2]) == VortexFilament(1.0,[v1,v2],[Segment(v1,v2),Segment(v2,v1)],[1,2],Int64[])
    @test VortexFilament(1.0,[v1,v2,v3]) == VortexFilament(1.0,[v1,v2,v3],[Segment(v1,v2),Segment(v2,v3),Segment(v3,v1)],[1,2,3],Int64[])
    @test VortexFilament(1.0,[v1,v2,v3],[2]) == VortexFilament(1.0,[v1,v2,v3],[Segment(v1,v2),Segment(v2,v3),Segment(v3,v1)],[1,3],[2])
end

@testset "isinf" begin
    v1 = [Inf,-1,0]
    v2 = [0,-1,0]
    v3 = [0,1,0]
    v4 = [Inf,1,0]
    vf = VortexFilament(1,[v1,v2,v3,v4])
    @test isinf(vf)

    v1 = [-Inf,1,0]
    v2 = [Inf,1,0]
    vf = VortexFilament(1,[v1,v2])
    @test isinf(vf)
end

@testset "issemiinf" begin
    v1 = rand(3)
    v2 = [Inf,rand(Float64),rand(Float64)]
    vf = VortexFilament(1,[v1,v2])
    @test issemiinf(vf)

    v1 = [rand(Float64),-Inf,rand(Float64)]
    v2 = rand(3)
    v3 = rand(3)
    vf = VortexFilament(1,[v1,v2,v3])
    @test issemiinf(vf)

    v1 = [Inf,-1,0]
    v2 = rand(3)
    v3 = rand(3)
    v4 = [Inf,1,0]
    vf = VortexFilament(1,[v1,v2,v3,v4])
    @test issemiinf(vf) == false
end

@testset "_closestpointonfinitesegmentaxis" begin
    v1 = [0,0,0]
    v2 = [1,0,0]
    @test VortexFilaments._closestpointonfinitesegmentaxis(Segment(v1,v2),[2,1,0]) == ([2,0,0],2)
    @test VortexFilaments._closestpointonfinitesegmentaxis(Segment(v1,v2),[-2,1,0]) == ([-2,0,0],-2)
end

@testset "infdir" begin
    v1 = rand(3)
    v2 = rand(3)
    v3 = rand(3)
    vtwoinf = [rand(Float64), Inf, Inf]
    vplusinfx = [Inf,rand(Float64),rand(Float64)]
    vplusinfy = [rand(Float64),Inf,rand(Float64)]
    vplusinfz = [rand(Float64),rand(Float64),Inf]
    vmininfx = [-Inf,rand(Float64),rand(Float64)]
    vmininfy = [rand(Float64),-Inf,rand(Float64)]
    vmininfz = [rand(Float64),rand(Float64),-Inf]
    @test_throws ErrorException VortexFilaments.infdir(VortexFilament(1.0,[v1,v2,vtwoinf]))
    @test_throws ErrorException VortexFilaments.infdir(VortexFilament(1.0,[vplusinfx,v2,vplusinfy]))
    @test VortexFilaments.infdir(VortexFilament(1.0,[v1,v2,vplusinfx])) == 1
    @test VortexFilaments.infdir(VortexFilament(1.0,[v1,v2,vplusinfy])) == 2
    @test VortexFilaments.infdir(VortexFilament(1.0,[v1,v2,vplusinfz])) == 3
    @test VortexFilaments.infdir(VortexFilament(1.0,[v1,v2,vmininfx])) == 1
    @test VortexFilaments.infdir(VortexFilament(1.0,[v1,v2,vmininfy])) == 2
    @test VortexFilaments.infdir(VortexFilament(1.0,[v1,v2,vmininfz])) == 3
    @test VortexFilaments.infdir(VortexFilament(1.0,[vplusinfx,v1,v2,vmininfx])) == 1
    @test VortexFilaments.infdir(VortexFilament(1.0,[vplusinfx,v1,v2,vplusinfx])) == 1
end

@testset "inducevelocity" begin
    v0 = [0,0,0]
    vplusinfx = [Inf,rand(Float64),rand(Float64)]
    vmininfx = [-Inf,rand(Float64),rand(Float64)]
    vplusinfx2 = [Inf,rand(Float64),rand(Float64)]
    vplusinfy = [rand(Float64),Inf,rand(Float64)]
    v1 = [-1,-1,0]
    v2 = [-1,1,0]
    v3 = [1,1,1]
    v4 = [1,-1,1]
    r = rand(Float64) + 1
    rvec = rand(3) .+ 1

    # Special cases for evaluation point
    s = Segment(v0,v3)
    @test inducevelocity(s,r*v0+(r-1)*v3) == [0.0,0.0,0.0] # evaluation point on segment
    @test inducevelocity(s,v0) == [0.0,0.0,0.0] # evaluation point on vertex 1
    @test inducevelocity(s,v3) == [0.0,0.0,0.0] # evaluation point on vertex 2
    @test inducevelocity(s,v0+(1+r)*v3) == [0.0,0.0,0.0] # evaluation point on segment axis but not on segment itself
    @test inducevelocity(Segment([0,1,2],[Inf,1,2]),[1,1,2]) == [0.0,0.0,0.0] # evaluation point on semi-infinite segment
    @test inducevelocity(Segment([0,1,2],[Inf,1,2]),[-1,1,2]) == [0.0,0.0,0.0] # evaluation point on semi-infinite segment axis but not on segment itself
    @test inducevelocity(Segment([-Inf,1,2],[Inf,1,2]),[1,1,2]) == [0.0,0.0,0.0] # evaluation point on infinite segment

    # Special segments
    s = Segment(v0,vplusinfx) # semi-infinite segment
    @test isapprox(norm(inducevelocity(s,[0.0,rvec[2],rvec[3]])), 1/(4π*norm([rvec[2],rvec[3]])), atol=1e-3)
    s = Segment(vmininfx,vplusinfx) # infinite segment
    centerpoint = 0.5*(vplusinfx + vmininfx)
    @test isapprox(norm(inducevelocity(s,rvec)), 1/(2π*norm([rvec[2]-centerpoint[2],rvec[3]-centerpoint[3]])), atol=1e-3)
    s = Segment(vmininfx,vplusinfy) # infinite segment with different infinite coordinates
    @test_throws ErrorException inducevelocity(s,rvec)
    s = Segment(vplusinfx,vplusinfx2) # infinite segment with non-opposite infinite vertices
    @test_throws ErrorException inducevelocity(s,rvec)

    # Horseshoe vortex
    b = rand(Float64)+1
    Γ = -rand(Float64)
    v1 = [Inf,-b/2,0]
    v2 = [0,-b/2,0]
    v3 = [0,b/2,0]
    v4 = [Inf,b/2,0]
    vf = VortexFilament(Γ,[v1,v2,v3,v4])
    downwash(Γ,b,y) = Γ/(4π)*b/((b/2)^2-y^2)
    yrange = range(-b/2,b/2,length=10)
    xevals = [[0.0,y,0.0] for y in yrange]
    @test isapprox(norm.(inducevelocity.(Ref(vf),xevals[2:end-1])), abs.(downwash.(Γ,b,yrange[2:end-1])), atol=1e-2)

    # Normal example
    v1 = [-1,-1,0]
    v2 = [-1,1,0]
    v3 = [1,1,1]
    v4 = [1,-1,1]

    vf = VortexFilament(1.0,[v1,v2,v3,v4])
    @test inducevelocity(vf,[-1,0,1]) ≈ [0.1566735365948611; 0.0; -0.12386284345966522]

    v1 = rand(3)
    v2 = rand(3)
    vtwodim = rand(2)
    @test_throws DimensionMismatch inducevelocity(VortexFilament(1.0,[v1,vtwodim,v2]),rand(3))

end
