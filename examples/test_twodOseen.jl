using FEMfunctions

import Gmsh: gmsh
using LinearAlgebra

using SparseArrays
using Test
using WriteVTK

include("makeAnnularMesh.jl")
include("twodOseen.jl")

μ = 1.0
ω =-1.0

x,eConn, innerNodes,innerX, outerNodes,outerX = makeAnnularMesh(1.0,2.0;lc=0.01)


#   Let's look at the mesh...
nNodes = size(x,2)
nElements = size(eConn,2)

xT = zeros(Float64,nNodes,2)   # xT = transpose(x[1:2,:])
eC = zeros(Int64,nElements,6)  # eC = transpose(eConn), converted to Int64
for i=1:nNodes
  xT[i,1] = x[1,i]
  xT[i,2] = x[2,i]
end
for i=1:nElements
  eC[i,1] = eConn[1,i]
  eC[i,2] = eConn[2,i]
  eC[i,3] = eConn[3,i]
  eC[i,4] = eConn[4,i]
  eC[i,5] = eConn[5,i]
  eC[i,6] = eConn[6,i]
end

advection = Array{Float64,2}(undef,nNodes,2)
for i=1:nNodes
  advection[i,1] = -ω*xT[i,2]
  advection[i,2] = ω*xT[i,1]
end

velocity,pressure,exactVelocity,exactPressure = twodOseen(xT,eC,innerNodes,outerNodes,advection,μ,ω)

#  To manually check the function saveFEMasVTK, enter the following line
#  and use ParaView:
# saveFEMasVTK("Oseen",xT,eC,["pressure"],pressure,["velocity"],velocity)
# saveFEMasVTK("OseenEx",xT,eC,["pressure"],exactPressure,["velocity"],exactVelocity)
# or combining these for convenience
# saveFEMasVTK("test.vtu",xT,eC,["pressure","exaPressure"],hcat(pressure,exactPressure),["velocity","exaVelocity"],hcat(velocity,exactVelocity))

M = twodMassMatrix(xT,eC)
errorU = velocity[:,1]-exactVelocity[:,1]
errorV = velocity[:,2]-exactVelocity[:,2]

pressure = pressure.+(exactPressure[1]-pressure[1])
errorP = pressure-exactPressure

@test sqrt( dot(errorU,(M*errorU)) + dot(errorV,(M*errorV)) ) < 1e-4
@test sqrt( dot(errorP,(M*errorP)) ) < 1e-2

