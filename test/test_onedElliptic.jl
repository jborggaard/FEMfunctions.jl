#  Solves a linear elliptic PDE in 1D 
  using FEMfunctions

  using LinearAlgebra
  using SparseArrays

  nElements = 40

  include("onedElliptic.jl")

  order = 1
  error = onedElliptic(nElements,order)
  @test error<1e-10

  order = 2
  error = onedElliptic(nElements,order)
  @test error<1e-10

  order = 3
  error = onedElliptic(nElements,order)
  @test error<1e-10
