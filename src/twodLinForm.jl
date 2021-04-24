function twodLinForm( Ff, test, wg )
#= """
  twodLinForm - routine to compute \int{ f*test }

  Author: Jeff Borggaard, Virginia Tech
  Version: 1.0

  Licensing:
     This code is distributed under the MIT license.
 
  Usage:
  ```julia
    F = twodLinForm( Ff, test, wg )
  ```

  Arguments:
  - `Ff`: Function values at the Gauss points
  - `test`: Matrix of test functions evaluated at the Gauss points 
            (dim: n_gauss, n_dof)
  - `wg`: Row vector of Gauss weights

  Output argument:
  - `F`: A numerical approximation to the integrals of phi * f over an element
""" =#

# n_dof = size(test,2)

# F = zeros(n_dof,1)
# for j=1:n_dof
#   F(j) = test(:,j)' * ( wg .* Ff )
# end
# F = sum(test'*(wg.*Ff),2)

  F = test'*(wg.*Ff)

  return F

end
