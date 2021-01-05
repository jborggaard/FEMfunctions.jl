using Test

@testset "FEMfunctions.jl" begin
  include("test_onedElliptic.jl")

  include("test_twodStokes.jl")

end
