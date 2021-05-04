using Random
using LinearAlgebra: norm

@testset "VortexFilament constructor" begin
    v1 = rand(3)
    v2 = rand(3)
    v3 = rand(3)
    vtwodim = rand(2)
    vtwoinf = [rand(Float64), Inf, Inf]
    vmininfx = [rand(Float64), rand(Float64), -Inf]
    vplusinfx = [rand(Float64), rand(Float64), Inf]
    vplusinfx2 = [rand(Float64), rand(Float64), Inf]
    vplusinfy = [rand(Float64), Inf, rand(Float64)]

    @test_throws AssertionError VortexFilament(1.0,[v1])
    @test_throws AssertionError VortexFilament(1.0,[v1,vtwodim,v3])
    @test_throws AssertionError VortexFilament(1.0,[v1,vtwoinf,v3])
    @test_throws AssertionError VortexFilament(1.0,[vplusinfx,vplusinfy])
    @test_throws AssertionError VortexFilament(1.0,[vplusinfx,vplusinfx2])
    @test VortexFilament(1.0,[vmininfx,vplusinfx]) == VortexFilament(1.0,[vmininfx,vplusinfx],[Segment(vmininfx,vplusinfx)],[1,2],Int64[])
    @test VortexFilament(1.0,[vplusinfx,vmininfx]) == VortexFilament(1.0,[vplusinfx,vmininfx],[Segment(vplusinfx,vmininfx)],[1,2],Int64[])
    @test VortexFilament(1.0,[v1,v2]) == VortexFilament(1.0,[v1,v2],[Segment(v1,v2),Segment(v2,v1)],[1,2],Int64[])
    @test VortexFilament(1.0,[v1,v2,v3]) == VortexFilament(1.0,[v1,v2,v3],[Segment(v1,v2),Segment(v2,v3),Segment(v3,v1)],[1,2,3],Int64[])
    @test VortexFilament(1.0,[v1,v2,v3],[2]) == VortexFilament(1.0,[v1,v2,v3],[Segment(v1,v2),Segment(v2,v3),Segment(v3,v1)],[1,3],[2])
    @test VortexFilament(1.0,[v1,v2,v3,v1]) == VortexFilament(1.0,[v1,v2,v3],[Segment(v1,v2),Segment(v2,v3),Segment(v3,v1)],[1,2,3],Int64[])

    r1 = rand(2)
    r2 = rand(2)
    for vinf1 in [[Inf,r1[1],r1[2]],[-Inf,r1[1],r1[2]],[r1[1],Inf,r1[2]],[r1[1],-Inf,r1[2]],[r1[1],r1[2],Inf],[r1[1],r1[2],-Inf]]
        @test VortexFilament(1.0,[v1,vinf1]) == VortexFilament(1.0,[v1,vinf1],[Segment(v1,vinf1)],[2],[1])
        @test VortexFilament(1.0,[vinf1,v1]) == VortexFilament(1.0,[vinf1,v1],[Segment(vinf1,v1)],[1],[2])
        @test VortexFilament(1.0,[v1,v2,vinf1]) == VortexFilament(1.0,[v1,v2,vinf1],[Segment(v1,v2),Segment(v2,vinf1)],[2,3],[1])
        @test VortexFilament(1.0,[vinf1,v1,v2]) == VortexFilament(1.0,[vinf1,v1,v2],[Segment(vinf1,v1),Segment(v1,v2)],[1,2],[3])
        @test_throws AssertionError VortexFilament(1.0,[v1,vinf1,v2])
        @test_throws AssertionError VortexFilament(1.0,[vinf1,v1,v2],[1])
        @test_throws AssertionError VortexFilament(1.0,[v1,v2,vinf1],[3])
        for vinf2 in [[Inf,r2[1],r2[2]],[-Inf,r2[1],r2[2]],[r2[1],Inf,r2[2]],[r2[1],-Inf,r2[2]],[r2[1],r2[2],Inf],[r2[1],r2[2],-Inf]]
            @test VortexFilament(1.0,[vinf1,v1,vinf2]) == VortexFilament(1.0,[vinf1,v1,vinf2],[Segment(vinf1,v1),Segment(v1,vinf2)],[1,2,3],Int64[])
            @test_throws AssertionError VortexFilament(1.0,[vinf1,v1,vinf2],[1])
            @test_throws AssertionError VortexFilament(1.0,[vinf1,v1,vinf2],[3])
            @test_throws AssertionError VortexFilament(1.0,[vinf1,v1,vinf2],[1,3])
        end
    end
end

@testset "_closestpointonsemiinfinitesegment" begin
    vfin = rand(3)
    pplus = vfin + rand(3)
    pmin = vfin - rand(3)
    for infdir in 1:3
        vinf = rand(3)
        vinf[infdir] = Inf

        s = Segment(vfin,vinf)
        closestpoint = vfin
        closestpoint[infdir] = pplus[infdir]
        @test VortexFilaments._closestpointonsemiinfinitesegment(s,pplus) == closestpoint
        @test VortexFilaments._closestpointonsemiinfinitesegment(s,pmin) == vfin

        s = Segment(vinf,vfin)
        closestpoint = vfin
        closestpoint[infdir] = pplus[infdir]
        @test VortexFilaments._closestpointonsemiinfinitesegment(s,pplus) == closestpoint
        @test VortexFilaments._closestpointonsemiinfinitesegment(s,pmin) == vfin

        s = Segment(vfin,-vinf)
        closestpoint = vfin
        closestpoint[infdir] = pmin[infdir]
        @test VortexFilaments._closestpointonsemiinfinitesegment(s,pmin) == closestpoint
        @test VortexFilaments._closestpointonsemiinfinitesegment(s,pplus) == vfin

        s = Segment(-vinf,vfin)
        closestpoint = vfin
        closestpoint[infdir] = pmin[infdir]
        @test VortexFilaments._closestpointonsemiinfinitesegment(s,pmin) == closestpoint
        @test VortexFilaments._closestpointonsemiinfinitesegment(s,pplus) == vfin
    end

end

@testset "_closestpointonfinitesegmentaxis" begin
    v1 = [0,0,0]
    v2 = [1,0,0]
    @test VortexFilaments._closestpointonfinitesegmentaxis(Segment(v1,v2),[2,1,0]) == ([2,0,0],2)
    @test VortexFilaments._closestpointonfinitesegmentaxis(Segment(v1,v2),[-2,1,0]) == ([-2,0,0],-2)
end

@testset "inducevelocity" begin
    v0 = [0,0,0]
    vplusinf = [Inf,0,0]
    vminusinf = [-Inf,0,0]
    v1 = [-1,-1,0]
    v2 = [-1,1,0]
    v3 = [1,1,1]
    v4 = [1,-1,1]
    r = rand(Float64) + 1
    r3 = rand(3) .+ 1

    # Special cases for evaluation point
    s = Segment(v0,v3)
    @test inducevelocity(s,r*v0+(r-1)*v3) == [0.0,0.0,0.0] # evaluation point on segment
    @test inducevelocity(s,v0) == [0.0,0.0,0.0] # evaluation point on vertex 1
    @test inducevelocity(s,v3) == [0.0,0.0,0.0] # evaluation point on vertex 2
    @test inducevelocity(s,v0+(1+r)*v3) == [0.0,0.0,0.0] # evaluation point on segment axis but not on segment itself
    s = Segment(v0,vplusinf) # semi-infinite segment
    @test inducevelocity(s,[1,0,0]) == [0.0,0.0,0.0] # evaluation point on segment
    @test inducevelocity(s,[-1,0,0]) == [0.0,0.0,0.0] # evaluation point on segment axis but not on segment itself
    s = Segment(vminusinf,vplusinf) # infinite segment
    @test inducevelocity(s,[1,0,0]) == [0.0,0.0,0.0] # evaluation point on segment

    # Special segments
    s = Segment(v0,vplusinf) # semi-infinite segment
    @test isapprox(norm(inducevelocity(s,[0.0,r3[2],r3[3]])), 1/(4π*norm([r3[2],r3[3]])), atol=1e-3)
    s = Segment(vminusinf,vplusinf) # infinite segment
    @test isapprox(norm(inducevelocity(s,r3)), 1/(2π*norm([r3[2],r3[3]])), atol=1e-3)

    # Horseshoe vortex
    b = rand(Float64)+1
    Γ = rand(Float64)
    v1 = [Inf,-b/2,0]
    v2 = [0,-b/2,0]
    v3 = [0,b/2,0]
    v4 = [Inf,b/2,0]
    vf = VortexFilament(Γ,[v1,v2,v3,v4])
    downwash(Γ,b,y) = Γ/(4π)*b/((b/2)^2-y^2)
    yrange = range(-b/2,b/2,length=10)
    xevals = [[0.0,y,0.0] for y in yrange]
    @test isapprox(norm.(inducevelocity.(Ref(vf),xevals[2:end-1])), downwash.(Γ,b,yrange[2:end-1]), atol=1e-3)

    # Normal example
    v1 = [-1,-1,0]
    v2 = [-1,1,0]
    v3 = [1,1,1]
    v4 = [1,-1,1]

    vf = VortexFilament(1.0,[v1,v2,v3,v4])
    @test inducevelocity(vf,[-1,0,1]) ≈ [0.1566735365948611; 0.0; -0.12386284345966522]



end
