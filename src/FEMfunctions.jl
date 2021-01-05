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

  export TriMesh_ProjectDerivatives
  export TriMesh_PromoteL2Q

  export twodBilinear
  export twodLinForm
  export twodMassMatrix
  export twodMesh
  export twodQuadratureRule
  export twodShape
  
  include("onedBilinear.jl")
  include("onedLinForm.jl")
  include("onedMesh.jl")
  include("onedQuadratureRule.jl")
  include("onedShape.jl")

  include("saveFEMasVTK.jl")

  include("TriMesh_ProjectDerivatives.jl")
  include("TriMesh_PromoteL2Q.jl")

  include("twodBilinear.jl")
  include("twodLinForm.jl")
  include("twodMassMatrix.jl")
  include("twodMesh.jl")
  include("twodQuadratureRule.jl")
  include("twodShape.jl")

end
