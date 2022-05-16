function makeAnnularMesh(rInner=1.0,rOuter=2.0;outDir=".",lc=0.1)
#  Uses Gmsh to create a triangular mesh in the annulus between two circles
#  defined by rInner and rOuter.  outDir is used to define the directory
#  containing the mesh file and lc is the nominal size of an element.
#
#  Usage:
#   x,eConn, nInner,xInner, nOuter,xOuter = 
#                                       makeAnnularMesh(rInner,rOuter;outDir,lc)
# 
#  An example triangular mesh used in FEMfunction examples
##

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 0) # change to 1 for Gmsh output

  gmsh.model.add("mixing")

  factory = gmsh.model.geo  # to shorten the gmsh.model.geo.* to factory.*
  # center point
  factory.addPoint(0.0,0.0,0.0, lc, 1)

  # four points along the inner circle
  factory.addPoint(rInner,0.0,0.0, lc, 2)
  factory.addPoint(0.0,rInner,0.0, lc, 3)
  factory.addPoint(-rInner,0.0,0.0, lc, 4)
  factory.addPoint(0.0,-rInner,0.0, lc, 5)

  # four points along the outer circle
  factory.addPoint(rOuter,0.0,0.0, lc, 6)
  factory.addPoint(0.0,rOuter,0.0, lc, 7)
  factory.addPoint(-rOuter,0.0,0.0, lc, 8)
  factory.addPoint(0.0,-rOuter,0.0, lc, 9)

  #  Assemble boundary points into curves following the right-hand-rule
  #  (traverse the boundary with your right hand toward the "outside")
  factory.addCircleArc(2,1,3, 11)        
  factory.addCircleArc(3,1,4, 12)
  factory.addCircleArc(4,1,5, 13)
  factory.addCircleArc(5,1,2, 14)

  factory.addCircleArc(6,1,7, 15)
  factory.addCircleArc(7,1,8, 16)
  factory.addCircleArc(8,1,9, 17)
  factory.addCircleArc(9,1,6, 18)

  factory.addCurveLoop([11,12,13,14],20)
  factory.addCurveLoop([15,16,17,18],30)
  factory.addPlaneSurface([-20,30],100)


  gmsh.model.addPhysicalGroup(1, [20], 40)  # 1D group is the inner circle  #40
  gmsh.model.addPhysicalGroup(1, [30], 50)  # 1D group is the outer circle, #50
  gmsh.model.setPhysicalName(1, 40, "Inner")
  gmsh.model.setPhysicalName(1, 50, "Outer")

  gmsh.model.addPhysicalGroup(2, [100], 1000)
#  gmsh.model.addPhysicalGroup(2, [2], 2)
  gmsh.model.setPhysicalName(2, 1000, "Mesh")

  factory.synchronize()

  gmsh.model.mesh.generate(2)
  gmsh.model.mesh.setOrder(2)  # quadratic elements

  ##nodeIds,nodeCoords,_ = gmsh.model.mesh.getNodes()
  ##x = reshape(nodeCoords,3,:)
  nodeTags1,xx1 = gmsh.model.mesh.getNodesForPhysicalGroup(2,1000)
  xx1 = reshape(xx1,3,:)
  nn = convert(Int64,maximum(nodeTags1))
  x = zeros(3,nn)
  for i=1:length(nodeTags1)
    x[:,nodeTags1[i]] = xx1[:,i]
  end

#  The following lines don't seem to work for Circular Arcs
#  nInner,xInner = gmsh.model.mesh.getNodesForPhysicalGroup(1,40) # 1D groups, #20 is the inner circle
#  nOuter,xOuter = gmsh.model.mesh.getNodesForPhysicalGroup(1,50) # 1D groups, #30 is the outer circle

#  Thus, we manually extract the nodes on the inner and outer circles.
  nInner11,xInner11 = gmsh.model.mesh.getNodes(1,11)
  nInner12,xInner12 = gmsh.model.mesh.getNodes(1,12)
  nInner13,xInner13 = gmsh.model.mesh.getNodes(1,13)
  nInner14,xInner14 = gmsh.model.mesh.getNodes(1,14)

  nOuter15,xOuter15 = gmsh.model.mesh.getNodes(1,15)
  nOuter16,xOuter16 = gmsh.model.mesh.getNodes(1,16)
  nOuter17,xOuter17 = gmsh.model.mesh.getNodes(1,17)
  nOuter18,xOuter18 = gmsh.model.mesh.getNodes(1,18)

  nInner = vcat(1:4,nInner11,nInner12,nInner13,nInner14)
  nOuter = vcat(5:8,nOuter15,nOuter16,nOuter17,nOuter18)

  sort!(nInner)
  sort!(nOuter)

  elemTypes1, elemTags1, elemConnectivity1 = gmsh.model.mesh.getElements(2,100)

  eC = reshape(elemConnectivity1[1],6,:)

  nElDOF,nElements = size(eC)
  eConn = Array{Int64,2}(undef,nElDOF,nElements)
  for i=1:nElements
    eConn[1,i] = eC[1,i]
    eConn[2,i] = eC[3,i]
    eConn[3,i] = eC[2,i]
    eConn[4,i] = eC[6,i]
    eConn[5,i] = eC[5,i]
    eConn[6,i] = eC[4,i]
  end

#  eConn = convert(Matrix{Int64},eConn)   # for debugging
#  nInner = convert(Matrix{Int64},nInner)
#  nOuter = convert(Matrix{Int64},nOuter)

  gmsh.write(outDir*"mixing.msh")
  gmsh.finalize()

  xInner = x[1:2,nInner]
  xOuter = x[1:2,nOuter]

 return x,eConn, nInner,xInner, nOuter,xOuter
end
