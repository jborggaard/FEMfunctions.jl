function twodOseen(x,eConn,innerNodes,outerNodes,advectionVelocity, 
                   μ::Float64=1.0, ω::Float64=1.0)
#
#  Solves the twod-dimensional steady-state Oseen equations
#       (v⋅∇)u + (u⋅∇)v = -∇p + ∇⋅(μ(∇u+∇u')) + f,
#                     0 = ∇⋅u
#  Ω is the domain between two concentric circles.  v is a prescribed 
#  advection velocity and f is used to generate a closed-form solution.
#
#  We prescribe zero Dirichlet conditions on the inner boundary,
#  and a uniform rotation of the outer boundary.
#
#  The parameter ϵ is used to implement a penalty method (computed using μ)
#
#  As a unit test, choose a value of f such that the solution u when
#  v = (-ωy,ωx) becomes u = (r-1)*(-ωy,ωx)/r and p = x^2+y^2.
#
#  v = (u,v) as the advection velocity, then u = (u,v) as the solution
#---------------------------------------------------------------------------78--

#  include("twodQuadratureRule.jl")
#  include("twodShape.jl")
#  include("twodBilinear.jl")
#  include("twodLinForm.jl")

  ϵ = 1e-8/μ   # penalty parameter
  rule  = 7    # points in quadrature formula

  function f(x::Array{Float64,2},μ,ω)
    d = size(x,1)
    ansX = zeros(Float64,d,1)
    ansY = zeros(Float64,d,1)
    for i=1:d
      rad = sqrt(x[i,1]^2+x[i,2]^2)
      ansX[i] = (2*x[i,1]*rad^5-3μ*ω*x[i,2]^3+ω^2*x[i,1]*(rad^2-rad^5) -
                  3*μ*ω*x[i,1]^2*x[i,2] + 4*μ*ω*x[i,2]*rad^2)/rad^5
      ansY[i] = (2*x[i,2]*rad^5+3μ*ω*x[i,1]^3+ω^2*x[i,2]*(rad^2-rad^5) +
                  3*μ*ω*x[i,1]*x[i,2]^2 - 4*μ*ω*x[i,1]*rad^2)/rad^5
    end
    return ansX,ansY
  end

  function uEx(x::Array{Float64,1},μ=1.0,ω=-1.0)
    rad = sqrt(x[1]^2+x[2]^2)
    factor = (rad-1.0)/rad

    ansX = -x[2]*ω*factor
    ansY =  x[1]*ω*factor
    ansP = rad^2
    
    return ansX, ansY, ansP
  end

  #  Get problem dimensions
  #-----------------------------------------------------------------------------
  nNodes = size(x,1)
  nElements = size(eConn,1)

#  innerNodes = sort!(innerNodes)  # assume they come in sorted
#  outerNodes = sort!(outerNodes)
  nInnerNodes = length(innerNodes)
  nOuterNodes = length(outerNodes)
#  nDirichlet = nInnerNodes+nOuterNodes   # use for Dirichlet BCs
  nDirichlet = nOuterNodes

  #  Set the index into equation numbers
  #-----------------------------------------------------------------------------
  ideU = zeros(Int64,nNodes,2)

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

  #  Integrate system matrices, element-by-element (each independently)
  #-------------------------------------------------------------------------78--
  r,s,w = twodQuadratureRule(rule)
  μ_g   = μ*ones(rule)
  ϵ_g   = ϵ*ones(rule)
  one   = ones(rule)

  nElDOF   = size(eConn,2)
  nElDOFsq = nElDOF*nElDOF
  nEntries = nElements*(2*nElDOF+3)*(2*nElDOF+3)

  II  = Array{Int64,1}(undef, nEntries)
  JJ  = Array{Int64,1}(undef, nEntries)
  AA  = Array{Float64,1}(undef, nEntries)
  b   = zeros(Float64,nUnk,1)

  for k=1:nElements
    nLocal = eConn[k,:][:]
    xLocal = x[nLocal,:]

    xg,wg,ϕ,ϕ_x,ϕ_y = twodShape( xLocal       , r, s, w )
    xg,wg,ψ,ψ_x,ψ_y = twodShape( xLocal[1:3,:], r, s, w )

    # evaluate forcing function at quadrature points
    fx_g, fy_g = f(xg,μ,ω)          

    u_g   = ϕ*advectionVelocity[nLocal,1]
    u_xg  = ϕ_x*advectionVelocity[nLocal,1]
    u_yg  = ϕ_y*advectionVelocity[nLocal,1]
    v_g   = ϕ*advectionVelocity[nLocal,2]
    v_xg  = ϕ_x*advectionVelocity[nLocal,2]
    v_yg  = ϕ_y*advectionVelocity[nLocal,2]

    A11Loc = -twodBilinear( 2*μ_g, ϕ_x, ϕ_x, wg) -
              twodBilinear(   μ_g, ϕ_y, ϕ_y, wg) -
              twodBilinear(   u_g, ϕ_x, ϕ  , wg) -
              twodBilinear(   v_g, ϕ_y, ϕ  , wg) -
              twodBilinear(  u_xg, ϕ  , ϕ  , wg)

    A12Loc = -twodBilinear(   μ_g, ϕ_x, ϕ_y, wg) -
              twodBilinear(  u_yg, ϕ  , ϕ  , wg)

    A21Loc = -twodBilinear(   μ_g, ϕ_y, ϕ_x, wg) -
              twodBilinear(  v_xg, ϕ  , ϕ  , wg)

    A22Loc = -twodBilinear(   μ_g, ϕ_x, ϕ_x, wg) -
              twodBilinear( 2*μ_g, ϕ_y, ϕ_y, wg) -
              twodBilinear(   u_g, ϕ_x, ϕ  , wg) -
              twodBilinear(   v_g, ϕ_y, ϕ  , wg) -
              twodBilinear(  v_yg, ϕ  , ϕ  , wg)
  
    B1Loc  =  twodBilinear(   one, ϕ_x, ψ  , wg)
    B2Loc  =  twodBilinear(   one, ϕ_y, ψ  , wg)

    MLoc   = -twodBilinear(   ϵ_g, ψ  , ψ  , wg)

    F1Loc  = twodLinForm( fx_g, ϕ, wg)
    F2Loc  = twodLinForm( fy_g, ϕ, wg)
    
    index = (k-1)*(2*nElDOF+3)^2  # compute base index (since k could be in any order)
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

  interiorNodes = convert(Array{Int64,1},1:nNodes)
  innerN = convert(Array{Int64,1},innerNodes)
  outerN = convert(Array{Int64,1},outerNodes)
  setdiff!(interiorNodes,innerN)
  setdiff!(interiorNodes,outerN)

  iB = vcat(innerN,outerN)
  nDirichlet = length(iB)
# knownIndexU = Array{Int64,1}(undef,nDirichlet)
# knownIndexV = Array{Int64,1}(undef,nDirichlet)
  knownIndexU = zeros(Int64,nDirichlet,1)
  knownIndexV = zeros(Int64,nDirichlet,1)

  dirichletU = zeros(Float64,nDirichlet,1)
  dirichletV = zeros(Float64,nDirichlet,1)
  for i=1:length(innerN)
    knownIndexU[i] = ideU[innerN[i],1] 
    knownIndexV[i] = ideU[innerN[i],2]
    uE,vE,pE = uEx(x[innerN[i],:],μ,ω)
    dirichletU[i] = uE
    dirichletV[i] = vE
  end

  ni = length(innerN)
  for i=1:length(outerN)
    idx = ni+i
    knownIndexU[idx] = ideU[outerN[i],1]
    knownIndexV[idx] = ideU[outerN[i],2]
    uE,vE,pE = uEx(x[outerN[i],:],μ,ω) 
    dirichletU[idx] = uE
    dirichletV[idx] = vE
  end

  nUnknowns = 2*length(interiorNodes) + nP
  unknownIndex = Array{UInt64,1}(undef,nUnknowns)

  global nUnknownCounter = 0
  for i ∈ interiorNodes
    global nUnknownCounter = nUnknownCounter+1
    unknownIndex[nUnknownCounter] = ideU[i,1]

    global nUnknownCounter = nUnknownCounter+1
    unknownIndex[nUnknownCounter] = ideU[i,2]
  end

  getP = zeros(Int64,3*nNodes,1)
  for nNode=1:nNodes
    if ideP[nNode]>0
      global nUnknownCounter = nUnknownCounter + 1
      unknownIndex[nUnknownCounter] = ideP[nNode]
      getP[ideP[nNode]] = nUnknownCounter
    end
  end

  rhs = -b[unknownIndex]-A[unknownIndex,knownIndexU[:]]*dirichletU-A[unknownIndex,knownIndexV[:]]*dirichletV

  uvp = A[unknownIndex,unknownIndex]\rhs

  velocity = zeros(Float64,nNodes,2)

  global nUnknownCounter = 0
  for i ∈ interiorNodes
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
  pressureExact = zeros(Float64,nNodes,1)
  for n=1:nNodes
    uE,vE,pE = uEx(x[n,:],μ,ω)
    velocityExact[n,1] = uE
    velocityExact[n,2] = vE
    pressureExact[n] = pE
  end

  # Compute the norm of the approximation error for the unit test
  M = twodMassMatrix(x,eConn)
#  vDiff = velocity-velocityExact
  C     = pressureExact[1]-pressure[1]
#  pDiff = pressure-pressureExact.+C
  pressureExact = pressureExact.-C
  return velocity,pressure,velocityExact,pressureExact
#  return sqrt(dot(vDiff[:,1],M,vDiff[:,1])) + 
#         sqrt(dot(vDiff[:,2],M,vDiff[:,2])) +
#         sqrt(dot(pDiff,M,pDiff))
end
