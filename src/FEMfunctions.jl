module FEMfunctions

  using LinearAlgebra
  using SparseArrays

  using WriteVTK

  export onedBilinear
  export onedLinForm
  export onedMesh
  export onedQuadratureRule
  export onedShape

  export saveFEMasVTK

  export TriMesh_ElementAdjacency
  export TriMesh_Interpolate
  export TriMesh_ProjectDerivatives
  export TriMesh_PromoteL2Q
  export TriMesh_Search
  export TriMesh_SquareMesh

  export twodBilinear
  export twodLinForm
  export twodMassMatrix
  export twodQuadratureRule
  export twodShape
  
  export TetMesh_CubeMesh
  export TetMesh_QuadratureRule

  export threedBilinear
  
  include("onedBilinear.jl")
  include("onedLinForm.jl")
  include("onedMesh.jl")
  include("onedQuadratureRule.jl")
  include("onedShape.jl")

  include("saveFEMasVTK.jl")

  include("TriMesh_ElementAdjacency.jl")
  include("TriMesh_Interpolate.jl")
  include("TriMesh_ProjectDerivatives.jl")
  include("TriMesh_PromoteL2Q.jl")
  include("TriMesh_Search.jl")
  include("TriMesh_SquareMesh.jl")

  include("twodBilinear.jl")
  include("twodLinForm.jl")
  include("twodMassMatrix.jl")
  include("twodQuadratureRule.jl")
  include("twodShape.jl")

  include("TetMesh_CubeMesh.jl")
  include("TetMesh_QuadratureRule.jl")

  include("threedBilinear.jl")

end
