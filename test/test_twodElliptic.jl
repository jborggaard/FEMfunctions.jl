using FEMfunctions

using LinearAlgebra
using SparseArrays

include("twodElliptic.jl")

x,eConn,u = twodElliptic(129)


function uExact(x)
  C = 0.2/π^2
  return C*sin(π*x[1])*sin(2*π*x[2])
end

nNodes = size(x,1)
uEx = zeros(Float64,nNodes)
for i=1:nNodes
  uEx[i] = uExact(x[i,:])
end

M = twodMassMatrix(x,eConn)

diff = u-uEx 
@test sqrt( dot(diff,(M*diff)) ) < 1e-8


