function twodAdvectionDiffusion(x,eConn,innerNodes,outerNodes,velocity, κ = 1.0,ω=-1.0)
#
#  Solves the advection-diffusion equation in 2D
#     - ∇⋅(κ∇θ) + v⋅∇θ = q,
#  Ω is an annular domain with outer radius 2 and inner radius 1.
#  We prescribe zero Dirichlet conditions to the outer boundary.
#
#  As a unit test, choose a value of q such that the solution θ when
#  v = (-ωy,ωx) becomes cos(2πr)-1 (which is zero at the outer boundary and has
#  zero flux on the inner boundary).
#---------------------------------------------------------------------------78--

#  include("twodQuadratureRule.jl")
#  include("twodShape.jl")
#  include("twodBilinear.jl")
#  include("twodLinForm.jl")

  rule  = 19;   # points in quadrature formula

  function q(x::Array{Float64,2},κ,ω)
    d = size(x,1)
    ans = zeros(Float64,d,1)
    for i=1:d
      rad = sqrt(x[i,1]^2+x[i,2]^2)
      ans[i] = 4*κ*π^2*cos(2*π*rad) + 2*κ*π*sin(2*π*rad)/rad
    end
    return ans
  end

  #  Get problem dimensions
  #-----------------------------------------------------------------------------
  nNodes = size(x,1)
  nElements = size(eConn,1)

  nInnerNodes = length(innerNodes)
  nOuterNodes = length(outerNodes)
#  nDirichlet = nInnerNodes+nOuterNodes   # use for Dirichlet BCs
  nDirichlet = nOuterNodes

  #  Set the index into equation numbers
  #-----------------------------------------------------------------------------
  ideθ = zeros(Int64,nNodes,1)

  global nUnk = 0
  for i=1:nNodes
    global nUnk = nUnk+1
    ideθ[i,1] = nUnk
  end

  #  Integrate system matrices, element-by-element
  #-----------------------------------------------------------------------------
  r,s,w = twodQuadratureRule(rule)
  κ_g   = κ*ones(rule)
  one   = ones(rule)

  nElDOF   = size(eConn,2)
  nElDOFsq = nElDOF*nElDOF
  nEntries = nElements*nElDOFsq
#  II  = SharedArray{Int64}(nEntries)
#  JJ  = SharedArray{Int64}(nEntries)
#  AA  = SharedArray{Float64}(nEntries)
#  b   = SharedArray{Float64}(nNodes)
  II  = Array{Int64,1}(undef, nEntries)
  JJ  = Array{Int64,1}(undef, nEntries)
  AA  = Array{Float64,1}(undef, nEntries)
  b   = zeros(Float64,nNodes,1)

#  @sync @distributed for k=1:nElements
  Threads.@threads for k=1:nElements
#  for k=1:nElements
    xg  = zeros(Float64,rule,2)
    wg  = zeros(Float64,rule)
    ϕ = zeros(Float64,rule,nElDOF)
    ϕ_x = zeros(Float64,rule,nElDOF)
    ϕ_y = zeros(Float64,rule,nElDOF)

    nLocal = eConn[k,:][:]
    xLocal = x[nLocal,:]

    xg,wg,ϕ,ϕ_x,ϕ_y = twodShape( xLocal, r, s, w )

    q_g   = q(xg,κ,ω)           # forcing function evaluated at quadrature points
    u_g   = ϕ*velocity[nLocal,1]
    v_g   = ϕ*velocity[nLocal,2]

    ALoc = twodBilinear( u_g, ϕ_x, ϕ  , wg) + twodBilinear( v_g, ϕ_y, ϕ  , wg) +
           twodBilinear( κ_g, ϕ_x, ϕ_x, wg) + twodBilinear( κ_g, ϕ_y, ϕ_y, wg) 
    FLoc = twodLinForm(  q_g, ϕ  ,      wg) 

    index = (k-1)*nElDOFsq  # compute base index (since k could be in any order)
    lDOF = ideθ[nLocal,1][:]

    for nt = 1:nElDOF
      nTest = ideθ[nLocal[nt],1]
      for nu = 1:nElDOF
        nUnkT = ideθ[nLocal[nu],1]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkT
        AA[index] = ALoc[nt,nu]
      end

      b[nTest] = b[nTest] + FLoc[nt]
    end
  end

  A = sparse(II,JJ,AA)

  interiorNodes = convert(Array{Int64,1},1:nNodes)
#  interiorNodes = setdiff(interiorNodes,innerNodes)
  outerN = convert(Array{Int64,1},outerNodes)
  setdiff!(interiorNodes,outerN)

#  knownIndexθ = [innerNodes;outerNodes]
  knownIndexθ = outerNodes

  dirichletθ = Array{Float64,1}(undef,nDirichlet)
  for i=1:nDirichlet
    dirichletθ[i] = 0.0
#    nO = knownIndexθ[i]
#    dirichletθ[i] = x[nO,1]
  end

  nUnknowns = length(interiorNodes)
  unKnownIndex = Array{UInt64,1}(undef,nUnknowns)

  global nUnknownCounter = 0
  for i ∈ interiorNodes
    global nUnknownCounter = nUnknownCounter+1
    unKnownIndex[nUnknownCounter] = i
  end

  θunk = A[unKnownIndex,unKnownIndex]\(b[unKnownIndex]-A[unKnownIndex,knownIndexθ[:]]*dirichletθ)

  θ = zeros(nNodes,1)

  global nUnknownCounter = 0
  for i ∈ interiorNodes
    global nUnknownCounter = nUnknownCounter+1
    θ[i] = θunk[nUnknownCounter]
  end

#  for i=1:nInnerNodes
#    θ[innerNodes[i]] = dirichletθ[i]
#  end

#  for i=1:nOuterNodes
#    θ[outerNodes[i]] = dirichletθ[nInnerNodes+i]
#  end
  for i=1:nOuterNodes
    θ[outerNodes[i]] = dirichletθ[i]
  end

  θexact = zeros(nNodes,1)
  for i=1:nNodes
    ri = sqrt(x[i,1]^2+x[i,2]^2)
    θexact[i] = cos(2*π*ri)-1.0
  end

  return θ,θexact
end
