function onedLinForm( Ff, test, wg )
#= """
  onedLinForm - routine to compute \int{ f*test }

  Author: Jeff Borggaard, Virginia Tech
  Version: 1.3

  Usage:
  ```julia
    F = onedLinForm( Ff, test, w_g )
  ```

  Arguments:
  - `Ff`: Function values at the Gauss points
  - `test`: Matrix of test functions evaluated at the Gauss points 
            (dim: n_gauss, n_dof)
  - `w_g`: Column vector of Gauss weights

  Output arguments:
  - `F`: the element's contribution to a right-hand-side
""" =#

  F = test'*(wg.*Ff)

  return F

end
