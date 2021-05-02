function threedBilinear( kernel, ϕ, test, w_g )
#= """
  threedBilinear.m - routine to compute \int{ kernel*ϕ*test }
                     (same as twodBilinear)
 
  Author: Jeff Borggaard, Virginia Tech
          part of FEMfunctions.jl
 
  Licensing:
     This code is distributed under the MIT license.

  Usage:
  ```julia
    M = threedBilinear(kernel, ϕ, test, w_g)
  ```
 
  Arguments:
  - `kernel`: Kernel function in the integral evaluated at the quadrature points
              (dim: nQuadrature, nDOF)
  - `ϕ`: Matrix of element test functions evaluated at the quadrature points 
              (dim: nQuadrature, nDOF)
  - `test`: Matrix of test functions evaluated at the quadrature points 
              (dim: nQuadrature, nDOF)
  - `wg`: Row vector of quadrature weights

  Output argument:
  - `M`: Numerical integral of kernel*ϕ*test
""" =#
#   [n_quadrature,n_row] = size(test)
#   [n_g1   ,n_col] = size(phi )
# 
#   M = zeros(n_row,n_col)
#   for i=1:n_row
#     for j=1:n_col
#        M(i,j) = ( w_g'    .* test(:,i)' ) * ( kernel .* phi(:,j) )
#     end
#   end

  wk = wg.*kernel

  M = test'*(LinearAlgebra.Diagonal(wk)*ϕ)

  return M
end
