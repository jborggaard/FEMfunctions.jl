module FEMfunctions

  include("onedBilinear.jl")
  include("onedLinForm.jl")
  include("onedMesh.jl")
  include("onedQuadratureRule.jl")
  include("onedShape.jl")

  include("twodBilinear.jl")
  include("twodLinForm.jl")
  include("twodMesh.jl")
  include("twodProjectDerivatives.jl")
  include("twodQuadratureRule.jl")
  include("twodShape.jl")

end
