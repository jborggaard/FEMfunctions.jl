function twodBilinear( kernel, ϕ, test, wg )
#= """
  twodBilinear - routine to compute \int{ kernel*ϕ*test }

  Author: Jeff Borggaard, Virginia Tech
          part of FEMfunctions.jl

  Licensing:
     This code is distributed under the MIT license.
 
  Usage:
  ```julia
    M = twodBilinear(kernel, ϕ, test, wg)
  ```

  Arguments:
  - `kernel`: Kernel function in the integral evaluated at the quadrature points
  - `ϕ`: Matrix of element test functions evaluated at the quadrature points 
           (dim: n_gauss, n_dof)
  - `test`: Matrix of test functions evaluated at the quadrature points 
            (dim: n_gauss, n_dof)
  - `wg`: Row vector of quadrature weights

  Output argument:
  - `M`: Numerical integral of kernel*phi*test
""" =#

#  M = test'*diagm( wg.*kernel )*ϕ


  wk = wg.*kernel

  M = test'*(LinearAlgebra.Diagonal(wk)*ϕ)

  return M

end
