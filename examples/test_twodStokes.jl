using FEMfunctions

using LinearAlgebra
using SparseArrays

include("twodStokes.jl")

@test twodStokes(25)<1e-10
