using FEMfunctions

import Gmsh: gmsh
using LinearAlgebra

using SparseArrays
#using Polynomials
using Printf
using WriteVTK

include("makeAnnularMesh.jl")
include("twodAdvectionDiffusion.jl")

κ = 1.0
ω = -1.0    # advection velocity parameter

x,eConn, innerNodes,innerX, outerNodes,outerX = makeAnnularMesh(1.0,2.0)

sort!(innerNodes);
sort!(outerNodes);


#   Let's look at the mesh...
nNodes = size(x,2)
nElements = size(eConn,2)

xT = zeros(Float64,nNodes,2)   # xT = transpose(x[1:2,:])
eC = zeros(Int64,nElements,6)  # eC = transpose(eConn), converted to Int64
for i=1:nNodes
  xT[i,1] = x[1,i]
  xT[i,2] = x[2,i]
end
for i=1:nElements#-nElementsSensor # the sensor elements are orientated correctly (why?)
  eC[i,1] = convert(Int64,eConn[1,i])
  eC[i,2] = convert(Int64,eConn[3,i])
  eC[i,3] = convert(Int64,eConn[2,i])
  eC[i,4] = convert(Int64,eConn[6,i])
  eC[i,5] = convert(Int64,eConn[5,i])
  eC[i,6] = convert(Int64,eConn[4,i])
end

velocity = Array{Float64,2}(undef,nNodes,2)
for i=1:nNodes
  velocity[i,1] = -ω*xT[i,2]
  velocity[i,2] = ω*xT[i,1]
end

temperature,exactTemperature = twodAdvectionDiffusion(xT,eC,innerNodes,outerNodes,velocity,κ,ω)

M = twodMassMatrix(xT,eC)
error = temperature-exactTemperature
#  The inner boundary isn't exactly a circle of radius 1, so we have a weaker error
@test sqrt( dot(error,(M*error)) ) < 2e-1

#  To manually check the function saveFEMasVTK, enter the following line
#  and use ParaView:
#  saveFEMasVTK("advDiff",xT,eC,["temperature"],temperature,["velocity"],velocity)
