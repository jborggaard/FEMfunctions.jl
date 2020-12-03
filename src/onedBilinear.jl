function onedBilinear( kernel, ϕ, test, w_g )
#function onedBilinear( kernel:: Array(Float64,1), 
#                       ϕ     :: Array(Float64,2), 
#                       test  :: Array(Float64,2),
#                       w_g   :: Array(Float64,1) )
#= """ 
  onedBilinear - routine to compute \int{ kernel*ϕ*test }
                 returns a matrix, 
                 rows correspond to different test functions, 
                 columns correspond to different unknowns (in ϕ)

  Author: Jeff Borggaard, Virginia Tech
  Version: 1.3

  Usage:
  ```julia
  M = onedBilinear(kernel, ϕ, test, w_g)
  ```

  Arguments:
  - `kernel`: Kernel function in the integral evaluated at the Gauss points
  - `ϕ`: Matrix of element test functions evaluated at the Gauss points 
           (dim: n_gauss, n_dof)
  - `test`: Matrix of test functions evaluated at the Gauss points 
            (dim: n_gauss, n_test)        
  - `w_g`: Column vector of Gauss weights

  Output arguments:
  - `M`: the current element's contribution to the system matrix
""" =#

  M = test'*diagm(kernel.*w_g)*ϕ

  return M

end
