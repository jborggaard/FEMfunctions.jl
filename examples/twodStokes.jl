function twodStokes(N=25)
#  Solves Stokes equation in 2D with Dirichlet boundary conditions
#     - ∇⋅(∇z+∇z') + ∇p = f,  Ω = (0,1)×(0,1),  homogeneous Dirichlet b.c.
#---------------------------------------------------------------------------78--

#  include("../src/twodMesh.jl")
#  include("../src/twodQuadratureRule.jl")
#  include("../src/twodShape.jl")
#  include("../src/twodBilinear.jl")
#  include("../src/twodLinForm.jl")
#  include("../src/twodMassMatrix.jl")
#  include("../src/TriMesh_PromoteL2Q.jl")

#  @everywhere include("twodShape.jl")
#  @everywhere include("twodBilinear.jl")
#  @everywhere include("twodLinForm.jl")

  #  Define problem parameters
#  @everywhere μ = 0.01
#  @everywhere ϵ = 0.0
  μ = 0.01
  ϵ = 0      # 1e-5/μ    # when the PDE solution is not in the FEM subspace 

  xMin = -1.0
  xMax =  1.0
  yMin = -1.0
  yMax =  1.0

  rule = 7   # points in quadrature formula


  #  Specify manufactured solutions for code verification
  #-----------------------------------------------------------------------------
  function uExact(x)
    return -x[2]*(1-x[2]^2)
  end

  function vExact(x)
    return  x[1]*(1-x[1]^2)
  end

  function pExact(x)
    return x[1]+x[2]
  end

#  @everywhere function fx(x::Array{Float64,2})
  function fx(x::Array{Float64,2})
    return 1.0-6.0*μ*x[1,2]
  end
#  @everywhere function fy(x::Array{Float64,2})
  function fy(x::Array{Float64,2})
    return 1.0+6.0*μ*x[1,1]
  end


  #  Create mesh then set equation numbers
  #-----------------------------------------------------------------------------
  nNodesX = N
  nNodesY = N
  x,eConn,iB = TriMesh_SquareMesh( xMin,xMax, 
                                   yMin,yMax, "quadratic", nNodesX,nNodesY )
  nNodes = size(x,1)
  nElements = size(eConn,1)


  #  Set the index into equation numbers
  #-----------------------------------------------------------------------------
  ideU = zeros(Int64,nNodes,2) 

  # easily streamlined... since boundaries are handled later
  global nUnk = 0
  for i=1:nNodes
    global nUnk = nUnk+1
    ideU[i,1] = nUnk

    global nUnk = nUnk+1
    ideU[i,2] = nUnk
  end

  ideP = zeros(Int64,nNodes)
  for n_el=1:nElements
    vertex = eConn[n_el,1:3]
    if ( ideP[vertex[1]]==0 )
      global nUnk = nUnk + 1
      ideP[vertex[1]] = nUnk
    end

    if ( ideP[vertex[2]]==0 )
      global nUnk = nUnk + 1
      ideP[vertex[2]] = nUnk
    end

    if ( ideP[vertex[3]]==0 )
      global nUnk = nUnk + 1
      ideP[vertex[3]] = nUnk
    end
  end
  nP = nUnk - 2*nNodes


  #  Integrate system matrices, element-by-element
  #-----------------------------------------------------------------------------
  r,s,w = twodQuadratureRule(rule)
  μ_g   = μ*ones(rule)
  ϵ_g   = ϵ*ones(rule)
  one   = ones(rule)

  nElDOF   = size(eConn,2)
  nElDOF2  = nElDOF*nElDOF
  nEntries = nElements*(2*nElDOF+3)*(2*nElDOF+3)
#  II  = SharedArray{Int64}(nEntries)
#  JJ  = SharedArray{Int64}(nEntries)
#  AA  = SharedArray{Float64}(nEntries)
#  b   = SharedArray{Float64}(2*nNodes)
  II  = Array{Int64,1}(undef, nEntries)
  JJ  = Array{Int64,1}(undef, nEntries)
  AA  = Array{Float64,1}(undef, nEntries)
  b   = zeros(Float64,nUnk,1)

#  @sync @distributed for k=1:nElements
  for k=1:nElements
#    xg  = zeros(Float64,rule,2)
#    wg  = zeros(Float64,rule)
#    phi = zeros(Float64,rule,nElDOF)
#    p_x = zeros(Float64,rule,nElDOF)
#    p_y = zeros(Float64,rule,nElDOF)

    nLocal = eConn[k,:][:]
    xLocal = x[nLocal,:]

    xg,wg,phi,phi_x,phi_y = twodShape( xLocal, r, s, w )
    xg,wg,psi,psi_x,psi_y = twodShape( xLocal[1:3,:], r, s, w )
    
    # evaluate forcing function at quadrature points
    fx_g = fx(xg)
    fy_g = fy(xg)

    A11Loc =-twodBilinear( 2*μ_g, phi_x, phi_x, wg) -
             twodBilinear(   μ_g, phi_y, phi_y, wg) 
    A12Loc =-twodBilinear(   μ_g, phi_x, phi_y, wg) 
    A21Loc =-twodBilinear(   μ_g, phi_y, phi_x, wg) 
    A22Loc =-twodBilinear(   μ_g, phi_x, phi_x, wg) -
             twodBilinear( 2*μ_g, phi_y, phi_y, wg) 
    B1Loc  = twodBilinear(   one, phi_x, psi  , wg) 
    B2Loc  = twodBilinear(   one, phi_y, psi  , wg) 
    MLoc   =-twodBilinear(   ϵ_g, psi  , psi  , wg) 
    F1Loc  = twodLinForm(   fx_g, phi  ,        wg) 
    F2Loc  = twodLinForm(   fy_g, phi  ,        wg) 

    index = (k-1)*(2*nElDOF+3)^2   # compute base index (since k could be in any order)
    lDOF = ideU[nLocal,1][:]
    # u-momentum equations
    for nt = 1:nElDOF
      nTest = ideU[nLocal[nt],1]
      for nu = 1:nElDOF
        nUnkU = ideU[nLocal[nu],1]
        nUnkV = ideU[nLocal[nu],2]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkU
        AA[index] = A11Loc[nt,nu]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkV
        AA[index] = A12Loc[nt,nu]
      end
      for np = 1:3  # special to this linear element
        nUnkP = ideP[nLocal[np]]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkP
        AA[index] = B1Loc[np,nt]
      end

      b[nTest] = b[nTest] + F1Loc[nt]
    end

    # v-momentum equations
    for nt = 1:nElDOF
      nTest = ideU[nLocal[nt],2]
      for nu = 1:nElDOF
        nUnkU = ideU[nLocal[nu],1]
        nUnkV = ideU[nLocal[nu],2]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkU
        AA[index] = A21Loc[nt,nu]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkV
        AA[index] = A22Loc[nt,nu]
      end
      for np = 1:3  # special to this linear element
        nUnkP = ideP[nLocal[np]]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkP
        AA[index] = B2Loc[np,nt]
      end

      b[nTest] = b[nTest] + F2Loc[nt]
    end

    # continuity equations
    for nt = 1:3
      nTest = ideP[nLocal[nt]]
      for nu = 1:nElDOF
        nUnkU = ideU[nLocal[nu],1]
        nUnkV = ideU[nLocal[nu],2]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkU
        AA[index] = B1Loc[nt,nu]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkV
        AA[index] = B2Loc[nt,nu]
      end
      for np = 1:3  # special to this linear element
        nUnkP = ideP[nLocal[np]]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkP
        AA[index] = MLoc[nt,np]
      end

    end
  end
  
  A = sparse(II,JJ,AA)


  #  Solve and apply the boundary conditions
  #-----------------------------------------------------------------------------
  nDirichlet = size(iB,1)

  sort!(iB)
  knownIndexU = [2*iB-ones(Int64,nDirichlet,1)]
  knownIndexV = [2*iB]

  dirichletU = zeros(Float64,nDirichlet,1)
  dirichletV = zeros(Float64,nDirichlet,1)
  for i = 1:length(iB)
    dirichletU[i] = uExact(x[iB[i],:])
    dirichletV[i] = vExact(x[iB[i],:])
  end   
  
  interiorNodes = 1:nNodes
  interiorNodes = setdiff(interiorNodes,iB)
#  sort!(interiorNodes)

  nUnknowns = 2*length(interiorNodes) + nP
  unknownIndex = Array{Int64,1}(undef,nUnknowns)
  uIndex = Array{Int64,1}(undef,2*length(interiorNodes))
  vIndex = Array{Int64,1}(undef,2*length(interiorNodes))
  uInverseIndex = Array{Int64,1}(undef,nNodes)
  vInverseIndex = Array{Int64,1}(undef,nNodes)

  global nUnknownCounter = 0
  for i ∈ interiorNodes
    global nUnknownCounter = nUnknownCounter+1
    unknownIndex[nUnknownCounter] = 2*i-1
    uIndex[nUnknownCounter] = i
    uInverseIndex[i]=nUnknownCounter

    global nUnknownCounter = nUnknownCounter+1
    unknownIndex[nUnknownCounter] = 2*i
    vIndex[nUnknownCounter] = i
    vInverseIndex[i]=nUnknownCounter
  end

  getP = zeros(Int64,3*nNodes,1)
  for nNode=1:nNodes
    if ideP[nNode]>0
      global nUnknownCounter = nUnknownCounter + 1
      unknownIndex[nUnknownCounter] = ideP[nNode]
      getP[ideP[nNode]] = nUnknownCounter
    end
  end

  rhs = -b[unknownIndex]-A[unknownIndex,knownIndexU[][:]]*dirichletU-A[unknownIndex,knownIndexV[][:]]*dirichletV

  uvp = A[unknownIndex,unknownIndex]\rhs

  velocity = zeros(Float64,nNodes,2)

  global nUnknownCounter = 0
  for i=interiorNodes
    global nUnknownCounter = nUnknownCounter+1
    velocity[i,1] = uvp[nUnknownCounter]

    global nUnknownCounter = nUnknownCounter + 1
    velocity[i,2] = uvp[nUnknownCounter]
  end

  for i=1:nDirichlet
    velocity[iB[i],1] = dirichletU[i]
    velocity[iB[i],2] = dirichletV[i]
  end

  pressure = zeros(Float64,nNodes,1)

  vv = eConn[:,1:3]
  vvv = vv[:]
  vertices = unique(sort(vvv))
  for i ∈ vertices
    pressure[i] = uvp[getP[ideP[i]]]
  end
  pressure = TriMesh_PromoteL2Q(pressure,eConn)


#  return x,eConn,velocity,pressure

  #  Compute the exact solution for the unit test
  velocityExact = zeros(Float64,nNodes,2)
  for n=1:nNodes
    velocityExact[n,1] = uExact(x[n,:])
    velocityExact[n,2] = vExact(x[n,:])
  end

  pressureExact = zeros(Float64,nNodes,1)
  for n=1:nNodes
    pressureExact[n] = pExact(x[n,:])
  end

  # Compute the norm of the approximation error for the unit test
  M = twodMassMatrix(x,eConn)
  vDiff = velocity-velocityExact
  C     = pressureExact[1]-pressure[1]
  pDiff = pressure-pressureExact.+C
  return sqrt(dot(vDiff[:,1],M,vDiff[:,1])) + 
         sqrt(dot(vDiff[:,2],M,vDiff[:,2])) +
         sqrt(dot(pDiff,M,pDiff))
end
