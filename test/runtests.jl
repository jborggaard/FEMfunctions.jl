#using FEMfunctions
using LinearAlgebra
using SparseArrays
using Test

@testset "FEMfunctions.jl" begin
  include("test_onedElliptic.jl")
  @test abs(test_onedElliptic())<1e-10

#  include("test_twodElliptic.jl")
#  @test abs(test_twodElliptic())<1e-10

  include("twodStokes.jl")
  @test abs(twodStokes(25))<0.01
end
