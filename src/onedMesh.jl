function onedMesh(xb, order::Int, numElem::Int)
#= """
  onedMesh    - Generate a one dimensional mesh with numElem elements.
                This routine returns elements of the same type as 
                the inputs xb, eConnb

  Author: Jeff Borggaard, Virginia Tech
          part of FEMfunctions.jl

  Licensing:
     This code is distributed under the MIT license.
 
  Usage:
  ```julia
    x, eConn, indexU, indexC = onedMesh(xb,order,numElem)
  ```

  Arguments:
  - `xb`: Endpoints of the interval
  - `order`: The order of the finite element mesh, 1-linear, 2-quadratic, etc.
  - `numElem`: Number of elements

  Output arguments:
  - `x`: Nodal coordinates of the mesh
  - `eConn`: Element connectivity of adapted mesh
  - `indexU`: Node numbers of unknowns
  - `indexC`: Node numbers of boundaries
""" =#

  #  Create the element connectivity
  eConn = zeros(Int64,numElem,order+1)
  for k=1:numElem
    if (order == 1)
      eConn[k,:] = [k, k+1]
    elseif (order == 2)
      eConn[k,:] = [2*k-1, 2*k, 2*k+1]
    elseif (order == 3)
      eConn[k,:] = [3*k-2, 3*k-1, 3*k, 3*k+1]
    else
      eConn[k,:] = order*k-(order-1):order*k+1
    end 
  end

  nNodes = eConn[end,end]
  x = collect(range(xb[1],xb[2],length=nNodes))

  indexU = [2:nNodes-1;]
  indexC = [1;nNodes]

  return x, eConn, indexU, indexC
end
