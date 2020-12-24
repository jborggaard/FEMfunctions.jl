module FEMfunctions

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
