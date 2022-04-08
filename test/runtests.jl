using Test
using Gmsh

@testset "FEMfunctions.jl" begin
  include("test_onedElliptic.jl")

  include("test_twodElliptic.jl")

  include("test_TriMesh.jl")

  include("../examples/test_twodAdvectionDiffusion.jl")

#  include("test_threedElliptic.jl")
#
#  include("test_TetMesh_QuadratureRule.jl")

end
