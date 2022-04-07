using FEMfunctions
using LinearAlgebra

x,eConn,iB = TriMesh_SquareMesh(0.0,1.0,0.0,1.0,"quadratic",5,5)

eAdjacency = TriMesh_ElementAdjacency(eConn)

xInterp = rand(20,2)

ff,eList = TriMesh_Interpolate(x,eConn,eAdjacency,x,xInterp)

@test norm(ff-xInterp)<1e-14
