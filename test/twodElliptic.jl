function twodElliptic(N=25, κ=1.0)
#
#  Solves a linear elliptic PDE on the unit square
#     - ∇⋅(κ∇u) = q
#---------------------------------------------------------------------------78--

#  include("twodMesh.jl")
#  include("twodQuadratureRule.jl")
#  @everywhere include("twodShape.jl")
#  @everywhere include("twodBilinear.jl")
#  @everywhere include("twodLinForm.jl")

  rule  = 7     # points in quadrature formula, adjust for different unit tests
                # options are 1, 3 (2nd degree), 7 (5th degree), 
                # 13 (7th degree), and 19 (9th degree).

  #  Set the verification example
  function uExact(x)
    C = 0.2/π^2
    return C*sin.(π*x[:,1])*sin.(2*π*x[:,2])
  end
#  f(x)=x[1]+x[2];            # right-hand side function definition.
#  @everywhere function f(x::Array{Float64,2})
  function f(x::Array{Float64,2})
    return sin.(π*x[:,1]).*sin.(2*π*x[:,2])
  end

  #  Specify the geometry/mesh
  #-----------------------------------------------------------------------------
  xMin = 0.0
  xMax = 1.0
  yMin = 0.0
  yMax = 1.0
  nNodesX = N 
  nNodesY = N
  x,eConn,iB = twodMesh( xMin,xMax, yMin,yMax, "quadratic", nNodesX,nNodesY )

  #  Get problem dimensions
  #-----------------------------------------------------------------------------
  nNodes = size(x,1)
  nElements = size(eConn,1)

  nDirichlet = length(iB)

  #  Set the index into equation numbers
  #-----------------------------------------------------------------------------
  ide = zeros(Int64,nNodes,1)
  
  global nUnk = 0
  for i=1:nNodes
    global nUnk = nUnk+1
    ide[i,1] = nUnk
  end

  # ide = 1:nNodes

  #  Integrate system matrices, element-by-element
  #-----------------------------------------------------------------------------
  r,s,w = twodQuadratureRule(rule)
  o     = ones(rule)
  κ_g   = κ*ones(rule)

  nElDOF   = size(eConn,2)
  nElDOF2  = nElDOF*nElDOF
  nEntries = nElements*nElDOF2
  II  = Array{Int64,1}(undef, nEntries)    #SharedArray(Int32,nEntries);
  JJ  = Array{Int64,1}(undef, nEntries)    #SharedArray(Int32,nEntries);
  AA  = Array{Float64,1}(undef, nEntries)  #SharedArray(Float64,nEntries);
  b   = zeros(Float64,nNodes,1)            #SharedArray(Float64,nNodes);

#  @sync @parallel for k=1:nElements
  for k=1:nElements
    #xg  = Array(Float64,rule,2);
    #wg  = Array(Float64,rule);
    #phi = Array(Float64,rule,nElDOF);
    #p_x = Array(Float64,rule,nElDOF);
    #p_y = Array(Float64,rule,nElDOF);
   
    nLocal = eConn[k,:][:]
    xLocal = x[nLocal,:]

    xg,wg,ϕ,ϕ_x,ϕ_y = twodShape( xLocal, r, s, w )
    fg   = f(xg)         # forcing function evaluated at quadrature points

    ALoc = twodBilinear( κ_g, ϕ_x, ϕ_x, wg ) + twodBilinear( κ_g, ϕ_y, ϕ_y, wg )
    bLoc = twodLinForm(   fg,      ϕ  , wg )

    index = (k-1)*nElDOF2    # compute base index (k could be in any order)
    lDOF = ide[nLocal,1][:]

    for nt = 1:nElDOF
      nTest = ide[nLocal[nt],1]
      for nu = 1:nElDOF
        nUnkU = ide[nLocal[nu],1]

        index = index + 1
        II[index] = nTest   # lDOF[nu]
        JJ[index] = nUnkU   # lDOF[nt]
        AA[index] = ALoc[nt,nu]
      end

      b[nTest] = b[nTest] + bLoc[nt] #b[lDOF[nt]] = b[lDOF[nt]] + bLoc[nt]
    end
  end

  A = sparse(II,JJ,AA)

  interiorNodes = 1:nNodes
  interiorNodes = setdiff(interiorNodes,iB)

  knownIndexU = iB

#  dirichletU = Array{Float64,1}(undef,nDirichlet)
#  for i=1:nDirichlet
#    dirichletU = uExact(x[i,:])
#  end

  u = zeros(Float64,nNodes)
  u[interiorNodes] = A[interiorNodes,interiorNodes]\b[interiorNodes]

  #print("$(u)")

  return x,eConn,u
end

