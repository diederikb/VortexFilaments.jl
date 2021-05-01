v1 = Vertex(-1,-1,0)
v2 = Vertex(-1,1,0)
v3 = Vertex(1,1,1)
v4 = Vertex(1,-1,1)
vertices = [v1,v2,v3,v4]

vl = VortexFilament(1.0,vertices)

xeval = [-1,0,1]

vel = inducevelocity(vl,xeval)

@test vel ≈ [0.1566735365948611; 0.0; -0.12386284345966522]
