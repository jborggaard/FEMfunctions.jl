using LinearAlgebra
using Test

include("../src/TetMesh_QuadratureRule.jl")

#  Tests a quadrature formula for tetrahedron.  Rules must produce the
#  stated accuracy for monomials when compared to the exact integral
#  over the unit tetrahedron:
#
#      int(int(int( monomial, z,0,1-x-y), y,0,1-x), x,0,1)
#

  # Test degree 0 monomials
  function pass0(r::Vector{Float64},s::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    err = zeros(Float64,1)
    err[1] = sum(w)-1.0/6.0
    return norm(err,Inf)<1e-15
  end

  # Test degree 1 monomials
  function pass1(r::Vector{Float64},s::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    err = zeros(Float64,3)
    err[1] = dot(r,w)-1.0/24.0
    err[2] = dot(s,w)-1.0/24.0
    err[3] = dot(t,w)-1.0/24.0
    return norm(err,Inf)<1e-15
  end

  # Test degree 2 monomials
  function pass2(r::Vector{Float64},s::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    err = zeros(Float64,6)
    err[1] = dot(r.*r,w)-1.0/60.0
    err[2] = dot(r.*s,w)-1.0/120.0
    err[3] = dot(r.*t,w)-1.0/120.0
    err[4] = dot(s.*s,w)-1.0/60.0
    err[5] = dot(s.*t,w)-1.0/120.0
    err[6] = dot(t.*t,w)-1.0/60.0
    return norm(err,Inf)<1e-15
  end

  # Test degree 3 monomials
  function pass3(r::Vector{Float64},s::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    err = zeros(Float64,10)
    err[ 1] = dot(r.*r.*r,w)-1.0/120.0
    err[ 2] = dot(r.*r.*s,w)-1.0/360.0
    err[ 3] = dot(r.*r.*t,w)-1.0/360.0
    err[ 4] = dot(r.*s.*s,w)-1.0/360.0
    err[ 5] = dot(r.*s.*t,w)-1.0/720.0
    err[ 6] = dot(r.*t.*t,w)-1.0/360.0
    err[ 7] = dot(s.*s.*s,w)-1.0/120.0
    err[ 8] = dot(s.*s.*t,w)-1.0/360.0
    err[ 9] = dot(s.*t.*t,w)-1.0/360.0
    err[10] = dot(t.*t.*t,w)-1.0/120.0
    return norm(err,Inf)<1e-15
  end

  # Test degree 4 monomials
  function pass4(r::Vector{Float64},s::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    err = zeros(Float64,15)
    err[ 1] = dot(r.*r.*r.*r,w)-1.0/210.0
    err[ 2] = dot(r.*r.*r.*s,w)-1.0/840.0
    err[ 3] = dot(r.*r.*r.*t,w)-1.0/840.0
    err[ 4] = dot(r.*r.*s.*s,w)-1.0/1260.0
    err[ 5] = dot(r.*r.*s.*t,w)-1.0/2520.0
    err[ 6] = dot(r.*r.*t.*t,w)-1.0/1260.0
    err[ 7] = dot(r.*s.*s.*s,w)-1.0/840.0
    err[ 8] = dot(r.*s.*s.*t,w)-1.0/2520.0
    err[ 9] = dot(r.*s.*t.*t,w)-1.0/2520.0
    err[10] = dot(r.*t.*t.*t,w)-1.0/840.0
    err[11] = dot(s.*s.*s.*s,w)-1.0/210.0
    err[12] = dot(s.*s.*s.*t,w)-1.0/840.0
    err[13] = dot(s.*s.*t.*t,w)-1.0/1260.0
    err[14] = dot(s.*t.*t.*t,w)-1.0/840.0
    err[15] = dot(t.*t.*t.*t,w)-1.0/210.0
    return norm(err,Inf)<1e-15
  end

  # Test degree 5 monomials
  function pass5(r::Vector{Float64},s::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    err = zeros(Float64,21)
    err[ 1] = dot(r.*r.*r.*r.*r,w)-1.0/336.0
    err[ 2] = dot(r.*r.*r.*r.*s,w)-1.0/1680.0
    err[ 3] = dot(r.*r.*r.*r.*t,w)-1.0/1680.0
    err[ 4] = dot(r.*r.*r.*s.*s,w)-1.0/3360.0
    err[ 5] = dot(r.*r.*r.*s.*t,w)-1.0/6720.0
    err[ 6] = dot(r.*r.*r.*t.*t,w)-1.0/3360.0
    err[ 7] = dot(r.*r.*s.*s.*s,w)-1.0/3360.0
    err[ 8] = dot(r.*r.*s.*s.*t,w)-1.0/10080.0
    err[ 9] = dot(r.*r.*s.*t.*t,w)-1.0/10080.0
    err[10] = dot(r.*r.*t.*t.*t,w)-1.0/3360.0
    err[11] = dot(r.*s.*s.*s.*s,w)-1.0/1680.0
    err[12] = dot(r.*s.*s.*s.*t,w)-1.0/6720.0
    err[13] = dot(r.*s.*s.*t.*t,w)-1.0/10080.0
    err[14] = dot(r.*s.*t.*t.*t,w)-1.0/6720.0
    err[15] = dot(r.*t.*t.*t.*t,w)-1.0/1680.0
    err[16] = dot(s.*s.*s.*s.*s,w)-1.0/336.0
    err[17] = dot(s.*s.*s.*s.*t,w)-1.0/1680.0
    err[18] = dot(s.*s.*s.*t.*t,w)-1.0/3360.0
    err[19] = dot(s.*s.*t.*t.*t,w)-1.0/3360.0
    err[20] = dot(s.*t.*t.*t.*t,w)-1.0/1680.0
    err[21] = dot(t.*t.*t.*t.*t,w)-1.0/336.0
    return norm(err,Inf)<1e-15
  end

  # Test degree 6 monomials
  function pass6(r::Vector{Float64},s::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    err = zeros(Float64,28)
    err[ 1] = dot(r.*r.*r.*r.*r.*r,w)-1.0/504.0
    err[ 2] = dot(r.*r.*r.*r.*r.*s,w)-1.0/3024.0
    err[ 3] = dot(r.*r.*r.*r.*r.*t,w)-1.0/3024.0
    err[ 4] = dot(r.*r.*r.*r.*s.*s,w)-1.0/7560.0
    err[ 5] = dot(r.*r.*r.*r.*s.*t,w)-1.0/15120.0
    err[ 6] = dot(r.*r.*r.*r.*t.*t,w)-1.0/7560.0
    err[ 7] = dot(r.*r.*r.*s.*s.*s,w)-1.0/10080.0
    err[ 8] = dot(r.*r.*r.*s.*s.*t,w)-1.0/30240.0
    err[ 9] = dot(r.*r.*r.*s.*t.*t,w)-1.0/30240.0
    err[10] = dot(r.*r.*r.*t.*t.*t,w)-1.0/10080.0
    err[11] = dot(r.*r.*s.*s.*s.*s,w)-1.0/7560.0
    err[12] = dot(r.*r.*s.*s.*s.*t,w)-1.0/30240.0
    err[13] = dot(r.*r.*s.*s.*t.*t,w)-1.0/45360.0
    err[14] = dot(r.*r.*s.*t.*t.*t,w)-1.0/30240.0
    err[15] = dot(r.*r.*t.*t.*t.*t,w)-1.0/7560.0
    err[16] = dot(r.*s.*s.*s.*s.*s,w)-1.0/3024.0
    err[17] = dot(r.*s.*s.*s.*s.*t,w)-1.0/15120.0
    err[18] = dot(r.*s.*s.*s.*t.*t,w)-1.0/30240.0
    err[19] = dot(r.*s.*s.*t.*t.*t,w)-1.0/30240.0
    err[20] = dot(r.*s.*t.*t.*t.*t,w)-1.0/15120.0
    err[21] = dot(r.*t.*t.*t.*t.*t,w)-1.0/3024.0
    err[22] = dot(s.*s.*s.*s.*s.*s,w)-1.0/504.0
    err[23] = dot(s.*s.*s.*s.*s.*t,w)-1.0/3024.0
    err[24] = dot(s.*s.*s.*s.*t.*t,w)-1.0/7560.0
    err[25] = dot(s.*s.*s.*t.*t.*t,w)-1.0/10080.0
    err[26] = dot(s.*s.*t.*t.*t.*t,w)-1.0/7560.0
    err[27] = dot(s.*t.*t.*t.*t.*t,w)-1.0/3024.0
    err[28] = dot(t.*t.*t.*t.*t.*t,w)-1.0/504.0
    return norm(err,Inf)<1e-15
  end

  # Test degree 7 monomials
  function pass7(r::Vector{Float64},s::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    err = zeros(Float64,36)
    err[ 1] = dot(r.*r.*r.*r.*r.*r.*r,w)-1.0/720.0
    err[ 2] = dot(r.*r.*r.*r.*r.*r.*s,w)-1.0/5040.0
    err[ 3] = dot(r.*r.*r.*r.*r.*r.*t,w)-1.0/5040.0
    err[ 4] = dot(r.*r.*r.*r.*r.*s.*s,w)-1.0/15120.0
    err[ 5] = dot(r.*r.*r.*r.*r.*s.*t,w)-1.0/30240.0
    err[ 6] = dot(r.*r.*r.*r.*r.*t.*t,w)-1.0/15120.0
    err[ 7] = dot(r.*r.*r.*r.*s.*s.*s,w)-1.0/25200.0
    err[ 8] = dot(r.*r.*r.*r.*s.*s.*t,w)-1.0/75600.0
    err[ 9] = dot(r.*r.*r.*r.*s.*t.*t,w)-1.0/75600.0
    err[10] = dot(r.*r.*r.*r.*t.*t.*t,w)-1.0/25200.0
    err[11] = dot(r.*r.*r.*s.*s.*s.*s,w)-1.0/25200.0
    err[12] = dot(r.*r.*r.*s.*s.*s.*t,w)-1.0/100800.0
    err[13] = dot(r.*r.*r.*s.*s.*t.*t,w)-1.0/151200.0
    err[14] = dot(r.*r.*r.*s.*t.*t.*t,w)-1.0/100800.0
    err[15] = dot(r.*r.*r.*t.*t.*t.*t,w)-1.0/25200.0
    err[16] = dot(r.*r.*s.*s.*s.*s.*s,w)-1.0/15120.0
    err[17] = dot(r.*r.*s.*s.*s.*s.*t,w)-1.0/75600.0
    err[18] = dot(r.*r.*s.*s.*s.*t.*t,w)-1.0/151200.0
    err[19] = dot(r.*r.*s.*s.*t.*t.*t,w)-1.0/151200.0
    err[20] = dot(r.*r.*s.*t.*t.*t.*t,w)-1.0/75600.0
    err[21] = dot(r.*r.*t.*t.*t.*t.*t,w)-1.0/15120.0
    err[22] = dot(r.*s.*s.*s.*s.*s.*s,w)-1.0/5040.0
    err[23] = dot(r.*s.*s.*s.*s.*s.*t,w)-1.0/30240.0
    err[24] = dot(r.*s.*s.*s.*s.*t.*t,w)-1.0/75600.0
    err[25] = dot(r.*s.*s.*s.*t.*t.*t,w)-1.0/100800.0
    err[26] = dot(r.*s.*s.*t.*t.*t.*t,w)-1.0/75600.0
    err[27] = dot(r.*s.*t.*t.*t.*t.*t,w)-1.0/30240.0
    err[28] = dot(r.*t.*t.*t.*t.*t.*t,w)-1.0/5040.0
    err[29] = dot(s.*s.*s.*s.*s.*s.*s,w)-1.0/720.0
    err[30] = dot(s.*s.*s.*s.*s.*s.*t,w)-1.0/5040.0
    err[31] = dot(s.*s.*s.*s.*s.*t.*t,w)-1.0/15120.0
    err[32] = dot(s.*s.*s.*s.*t.*t.*t,w)-1.0/25200.0
    err[33] = dot(s.*s.*s.*t.*t.*t.*t,w)-1.0/25200.0
    err[34] = dot(s.*s.*t.*t.*t.*t.*t,w)-1.0/15120.0
    err[35] = dot(s.*t.*t.*t.*t.*t.*t,w)-1.0/5040.0
    err[36] = dot(t.*t.*t.*t.*t.*t.*t,w)-1.0/720.0
    return norm(err,Inf)<1e-15
  end

  # Test degree 8 monomials
  function pass8(r::Vector{Float64},s::Vector{Float64},t::Vector{Float64},w::Vector{Float64})
    err = zeros(Float64,45)
    err[ 1] = dot(r.*r.*r.*r.*r.*r.*r.*r,w)-1.0/990.0
    err[ 2] = dot(r.*r.*r.*r.*r.*r.*r.*s,w)-1.0/7920.0
    err[ 3] = dot(r.*r.*r.*r.*r.*r.*r.*t,w)-1.0/7920.0
    err[ 4] = dot(r.*r.*r.*r.*r.*r.*s.*s,w)-1.0/27720.0
    err[ 5] = dot(r.*r.*r.*r.*r.*r.*s.*t,w)-1.0/55440.0
    err[ 6] = dot(r.*r.*r.*r.*r.*r.*t.*t,w)-1.0/27720.0
    err[ 7] = dot(r.*r.*r.*r.*r.*s.*s.*s,w)-1.0/55440.0
    err[ 8] = dot(r.*r.*r.*r.*r.*s.*s.*t,w)-1.0/166320.0
    err[ 9] = dot(r.*r.*r.*r.*r.*s.*t.*t,w)-1.0/166320.0
    err[10] = dot(r.*r.*r.*r.*r.*t.*t.*t,w)-1.0/55440.0
    err[11] = dot(r.*r.*r.*r.*s.*s.*s.*s,w)-1.0/69300.0
    err[12] = dot(r.*r.*r.*r.*s.*s.*s.*t,w)-1.0/277200.0
    err[13] = dot(r.*r.*r.*r.*s.*s.*t.*t,w)-1.0/415800.0
    err[14] = dot(r.*r.*r.*r.*s.*t.*t.*t,w)-1.0/277200.0
    err[15] = dot(r.*r.*r.*r.*t.*t.*t.*t,w)-1.0/69300.0
    err[16] = dot(r.*r.*r.*s.*s.*s.*s.*s,w)-1.0/55440.0
    err[17] = dot(r.*r.*r.*s.*s.*s.*s.*t,w)-1.0/277200.0
    err[18] = dot(r.*r.*r.*s.*s.*s.*t.*t,w)-1.0/554400.0
    err[19] = dot(r.*r.*r.*s.*s.*t.*t.*t,w)-1.0/554400.0
    err[20] = dot(r.*r.*r.*s.*t.*t.*t.*t,w)-1.0/277200.0
    err[21] = dot(r.*r.*r.*t.*t.*t.*t.*t,w)-1.0/55440.0
    err[22] = dot(r.*r.*s.*s.*s.*s.*s.*s,w)-1.0/27720.0
    err[23] = dot(r.*r.*s.*s.*s.*s.*s.*t,w)-1.0/166320.0
    err[24] = dot(r.*r.*s.*s.*s.*s.*t.*t,w)-1.0/415800.0
    err[25] = dot(r.*r.*s.*s.*s.*t.*t.*t,w)-1.0/554400.0
    err[26] = dot(r.*r.*s.*s.*t.*t.*t.*t,w)-1.0/415800.0
    err[27] = dot(r.*r.*s.*t.*t.*t.*t.*t,w)-1.0/166320.0
    err[28] = dot(r.*r.*t.*t.*t.*t.*t.*t,w)-1.0/27720.0
    err[29] = dot(r.*s.*s.*s.*s.*s.*s.*s,w)-1.0/7920.0
    err[30] = dot(r.*s.*s.*s.*s.*s.*s.*t,w)-1.0/55440.0
    err[31] = dot(r.*s.*s.*s.*s.*s.*t.*t,w)-1.0/166320.0
    err[32] = dot(r.*s.*s.*s.*s.*t.*t.*t,w)-1.0/277200.0
    err[33] = dot(r.*s.*s.*s.*t.*t.*t.*t,w)-1.0/277200.0
    err[34] = dot(r.*s.*s.*t.*t.*t.*t.*t,w)-1.0/166320.0
    err[35] = dot(r.*s.*t.*t.*t.*t.*t.*t,w)-1.0/55440.0
    err[36] = dot(r.*t.*t.*t.*t.*t.*t.*t,w)-1.0/7920.0
    err[37] = dot(s.*s.*s.*s.*s.*s.*s.*s,w)-1.0/990.0
    err[38] = dot(s.*s.*s.*s.*s.*s.*s.*t,w)-1.0/7920.0
    err[39] = dot(s.*s.*s.*s.*s.*s.*t.*t,w)-1.0/27720.0
    err[40] = dot(s.*s.*s.*s.*s.*t.*t.*t,w)-1.0/55440.0
    err[41] = dot(s.*s.*s.*s.*t.*t.*t.*t,w)-1.0/69300.0
    err[42] = dot(s.*s.*s.*t.*t.*t.*t.*t,w)-1.0/55440.0
    err[43] = dot(s.*s.*t.*t.*t.*t.*t.*t,w)-1.0/27720.0
    err[44] = dot(s.*t.*t.*t.*t.*t.*t.*t,w)-1.0/7920.0
    err[45] = dot(t.*t.*t.*t.*t.*t.*t.*t,w)-1.0/990.0
    return norm(err,Inf)<1e-15
  end




  r,s,t,w = TetMesh_QuadratureRule(1)
  @test pass0(r,s,t,w)
  @test pass1(r,s,t,w)

#  r,s,t,w = TetMesh_QuadratureRule(4)
#  @test pass0(r,s,t,w)
#  @test pass1(r,s,t,w)
#  @test pass2(r,s,t,w)

  r,s,t,w = TetMesh_QuadratureRule(5)
  @test pass0(r,s,t,w)
  @test pass1(r,s,t,w)
  @test pass2(r,s,t,w)
  @test pass3(r,s,t,w)

  r,s,t,w = TetMesh_QuadratureRule(11)
  @test pass0(r,s,t,w)
  @test pass1(r,s,t,w)
  @test pass2(r,s,t,w)
  @test pass3(r,s,t,w)
  @test pass4(r,s,t,w)

  r,s,t,w = TetMesh_QuadratureRule(15)
  @test pass0(r,s,t,w)
  @test pass1(r,s,t,w)
  @test pass2(r,s,t,w)
  @test pass3(r,s,t,w)
  @test pass4(r,s,t,w)
  @test pass5(r,s,t,w)

  r,s,t,w = TetMesh_QuadratureRule(24)
  @test pass0(r,s,t,w)
  @test pass1(r,s,t,w)
  @test pass2(r,s,t,w)
  @test pass3(r,s,t,w)
  @test pass4(r,s,t,w)
  @test pass5(r,s,t,w)
  @test pass6(r,s,t,w)

  r,s,t,w = TetMesh_QuadratureRule(31)
  @test pass0(r,s,t,w)
  @test pass1(r,s,t,w)
  @test pass2(r,s,t,w)
  @test pass3(r,s,t,w)
  @test pass4(r,s,t,w)
  @test pass5(r,s,t,w)
  @test pass6(r,s,t,w)
  @test pass7(r,s,t,w)

  r,s,t,w = TetMesh_QuadratureRule(45)
  @test pass0(r,s,t,w)
  @test pass1(r,s,t,w)
  @test pass2(r,s,t,w)
  @test pass3(r,s,t,w)
  @test pass4(r,s,t,w)
  @test pass5(r,s,t,w)
  @test pass6(r,s,t,w)
  @test pass7(r,s,t,w)
  @test pass8(r,s,t,w)
