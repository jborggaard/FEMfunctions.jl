function TriMesh_PromoteL2Q(p,eConn)
#= 
  TriMesh_PromoteL2Q
   Interpolate values at the midside nodes from the vertices in a
   triangular mesh.  Used to promote a linear to a quadratic representation
   of a finite element function.
  
  Author: Jeff Borggaard, Virginia Tech
          part of FEMfunctions.jl

  Licensing:
     This code is distributed under the MIT license.
=# 

  nElements = size(eConn,1)

  for nEl=1:nElements
    localNodes = eConn[nEl,:]
    p[ localNodes[4] ] = ( p[localNodes[1]] + p[localNodes[2]] )/2.0
    p[ localNodes[5] ] = ( p[localNodes[2]] + p[localNodes[3]] )/2.0
    p[ localNodes[6] ] = ( p[localNodes[3]] + p[localNodes[1]] )/2.0
  end

  return p
end
