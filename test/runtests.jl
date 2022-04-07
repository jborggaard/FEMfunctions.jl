using Test
using Gmsh

@testset "FEMfunctions.jl" begin
  include("test_onedElliptic.jl")

  include("test_twodElliptic.jl")

  include("test_TriMesh.jl")

  include("../examples/test_twodAdvectionDiffusion.jl")
end
