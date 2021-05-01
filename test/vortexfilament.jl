using Random

@testset "VortexFilament constructor" begin
    v1 = Vertex(rand(3))
    v2 = Vertex(rand(3))
    v3 = Vertex(rand(3))

    @test_throws AssertionError VortexFilament(1.0,[v1])
    @test VortexFilament(1.0,[v1,v2]) == VortexFilament(1.0,[v1,v2],[Segment(v1,v2),Segment(v2,v1)],[1,2],Int64[])
    @test VortexFilament(1.0,[v1,v2,v3]) == VortexFilament(1.0,[v1,v2,v3],[Segment(v1,v2),Segment(v2,v3),Segment(v3,v1)],[1,2,3],Int64[])
    @test VortexFilament(1.0,[v1,v2,v3],[2]) == VortexFilament(1.0,[v1,v2,v3],[Segment(v1,v2),Segment(v2,v3),Segment(v3,v1)],[1,3],[2])
    @test VortexFilament(1.0,[v1,v2,v3,v1]) == VortexFilament(1.0,[v1,v2,v3],[Segment(v1,v2),Segment(v2,v3),Segment(v3,v1)],[1,2,3],Int64[])

    # r1 = rand(2)
    # r2 = rand(2)
    r1 = [1,2]
    r2 = [3,4]
    for vinf1 in [Vertex(Inf,r1[1],r1[2]),Vertex(-Inf,r1[1],r1[2]),Vertex(r1[1],Inf,r1[2]),Vertex(r1[1],-Inf,r1[2]),Vertex(r1[1],r1[2],Inf),Vertex(r1[1],r1[2],-Inf)]
        @test VortexFilament(1.0,[v1,vinf1]) == VortexFilament(1.0,[v1,vinf1],[Segment(v1,vinf1)],[2],[1])
        @test VortexFilament(1.0,[vinf1,v1]) == VortexFilament(1.0,[vinf1,v1],[Segment(vinf1,v1)],[1],[2])
        @test VortexFilament(1.0,[v1,v2,vinf1]) == VortexFilament(1.0,[v1,v2,vinf1],[Segment(v1,v2),Segment(v2,vinf1)],[2,3],[1])
        @test VortexFilament(1.0,[vinf1,v1,v2]) == VortexFilament(1.0,[vinf1,v1,v2],[Segment(vinf1,v1),Segment(v1,v2)],[1,2],[3])
        @test_throws AssertionError VortexFilament(1.0,[v1,vinf1,v2])
        @test_throws AssertionError VortexFilament(1.0,[vinf1,v1,v2],[1])
        @test_throws AssertionError VortexFilament(1.0,[v1,v2,vinf1],[3])
        for vinf2 in [Vertex(Inf,r2[1],r2[2]),Vertex(-Inf,r2[1],r2[2]),Vertex(r2[1],Inf,r2[2]),Vertex(r2[1],-Inf,r2[2]),Vertex(r2[1],r2[2],Inf),Vertex(r2[1],r2[2],-Inf)]
            @test VortexFilament(1.0,[vinf1,vinf2]) == VortexFilament(1.0,[vinf1,vinf2],[Segment(vinf1,vinf2)],[1,2],Int64[])
            @test VortexFilament(1.0,[vinf1,v1,vinf2]) == VortexFilament(1.0,[vinf1,v1,vinf2],[Segment(vinf1,v1),Segment(v1,vinf2)],[1,2,3],Int64[])
            @test_throws AssertionError VortexFilament(1.0,[vinf1,v1,vinf2],[1])
            @test_throws AssertionError VortexFilament(1.0,[vinf1,v1,vinf2],[3])
            @test_throws AssertionError VortexFilament(1.0,[vinf1,v1,vinf2],[1,3])
        end
    end
end

@testset "inducevelocity" begin
    v1 = Vertex(-1,-1,0)
    v2 = Vertex(-1,1,0)
    v3 = Vertex(1,1,1)
    v4 = Vertex(1,-1,1)
    vertices = [v1,v2,v3,v4]

    vl = VortexFilament(1.0,vertices)

    xeval = [-1,0,1]

    vel = inducevelocity(vl,xeval)

    @test vel ≈ [0.1566735365948611; 0.0; -0.12386284345966522]
end
