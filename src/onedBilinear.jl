function onedBilinear( kernel, ϕ, test, wg )
#function onedBilinear( kernel:: Array(Float64,1), 
#                       ϕ     :: Array(Float64,2), 
#                       test  :: Array(Float64,2),
#                       wg    :: Array(Float64,1) )
#= """ 
  onedBilinear - routine to compute (kernel*ϕ,test) = \int{ kernel*ϕ*test }
                 returns a matrix, 
                 rows correspond to different test functions, 
                 columns correspond to different local unknowns (in ϕ)

  Author: Jeff Borggaard, Virginia Tech
          part of FEMfunctions.jl

  Licensing:
     This code is distributed under the MIT license.
 
  Usage:
  ```julia
    M = onedBilinear(kernel, ϕ, test, wg)
  ```

  Arguments:
  - `kernel`: Kernel function in the integral evaluated at the Gauss points
  - `ϕ`: Matrix of element test functions evaluated at the Gauss points 
           (dim: n_gauss, n_dof)
  - `test`: Matrix of test functions evaluated at the Gauss points 
            (dim: n_gauss, n_test)        
  - `wg`: Column vector of Gauss weights

  Output arguments:
  - `M`: the current element's contribution to the system matrix
""" =#

  wk = wg.*kernel

  M = test'*(LinearAlgebra.Diagonal(wk)*ϕ)
 
  return M

end
