function makeAnnularMesh(rinner=1.0,router=2.0;outDir=".",lc=0.01)

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 0) # change to 1 for Gmsh output

  gmsh.model.add("mixing")

  n=100
  global theta = range(0.0,stop=2.0*pi,length=n+1)

  for i=1:n
    gmsh.model.geo.addPoint(rinner*cos(theta[i]),rinner*sin(theta[i]),
                            0.0, lc, i)
  end

  for i=1:n
    gmsh.model.geo.addPoint(router*cos(theta[i]),router*sin(theta[i]),
                            0.0, 2*lc,n+i)
  end

  #  Assemble boundary points into curves following the right-hand-rule
  #  (traverse the boundary with your right hand toward the "outside")
  gmsh.model.geo.addBSpline([1;collect(n:-1:1)],1)        
  gmsh.model.geo.addSpline([collect((n+1):(2*n));n+1],2)

  gmsh.model.geo.addCurveLoop([1],1)
  gmsh.model.geo.addCurveLoop([2],2)
  gmsh.model.geo.addPlaneSurface([1,2],1)

  gmsh.model.addPhysicalGroup(1, [1], 1)
  gmsh.model.addPhysicalGroup(1, [2], 2)
  gmsh.model.setPhysicalName(1, 1, "Inner")
  gmsh.model.setPhysicalName(1, 2, "Outer")

  gmsh.model.addPhysicalGroup(2, [1], 1)
  gmsh.model.addPhysicalGroup(2, [2], 2)
  gmsh.model.setPhysicalName(2, 1, "Mesh")

  gmsh.model.geo.synchronize()

  gmsh.model.mesh.generate(2)
  gmsh.model.mesh.setOrder(2)  # quadratic elements

  ##nodeIds,nodeCoords,_ = gmsh.model.mesh.getNodes()
  ##x = reshape(nodeCoords,3,:)
  nodeTags1,xx1 = gmsh.model.mesh.getNodesForPhysicalGroup(2,1)
  xx1 = reshape(xx1,3,:)
  nn = convert(Int64,maximum(nodeTags1))
  x = zeros(3,nn)
  for i=1:length(nodeTags1)
    x[:,nodeTags1[i]] = xx1[:,i]
  end

  nInner,xInner = gmsh.model.mesh.getNodesForPhysicalGroup(1,1) # 1D groups, #1 is the inner circle
  nOuter,xOuter = gmsh.model.mesh.getNodesForPhysicalGroup(1,2) # 1D groups, #2 is the outer circle

  elemTypes1, elemTags1, elemConnectivity1 = gmsh.model.mesh.getElements(2,1)

  eConn = reshape(elemConnectivity1[1],6,:)

#  eC = convert(Array{Int64,2},eConn)   # for debugging

#  gmsh.write(outDir*"mixing.msh")
  gmsh.finalize()

 return x,eConn, nInner,xInner, nOuter,xOuter
end
