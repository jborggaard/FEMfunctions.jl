using Test

@testset "FEMfunctions.jl" begin
  include("test_onedElliptic.jl")

  include("test_twodElliptic.jl")

  include("test_TriMesh.jl")

  include("test_twodAdvectionDiffusion.jl")
end
