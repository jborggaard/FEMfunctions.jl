function onedShape(x,r,w)
#= """
  oned_shape.m - computes test functions and derivatives for a Lagrange
                 C0 element given element coordinates and Gauss points.
                 (assumes all nodes are uniformly distributed in the
                 element)

  Author: Jeff Borggaard, Virginia Tech
  Version: 1.3

  Usage:
  ```julia
    x_g, w_g, ϕ, ϕ_x, ϕ_xx = onedShape(x,r,w)
  ```

  Arguments:
  - `x`: Coordinates of the element nodes
  - `r`: Coordinates of Gauss points in (-1,1)
  - `w`: Gauss weights associated with r

  Output arguments:
  - `x_g`: Coordinates of Gauss points in the element
  - `w_g`: Gauss weights scaled by the element Jacobian
  - `ϕ`: Value of element shape functions at x_g
  - `ϕ_x`: First spatial derivatives of ϕ
  - `ϕ_xx`: Second spatial derivatives of ϕ
""" =#

  n  = size(x,1)
  ng = size(r,1)

  xg  = zeros(Float64,ng,n)
  ϕ   = zeros(Float64,ng,n)
  ϕ_x = zeros(Float64,ng,n)
  ϕ_xx = zeros(Float64,ng,n)

  if n == 2            
    # Transform coordinates for linear elements
    c0 = ( x[n]-x[1] )/2.0
    c1 = ( x[n]+x[1] )/2.0

    x_g = c0.*r .+ c1

    ϕ[:,1] = ( 1-r )/2.0
    ϕ[:,2] = ( 1+r )/2.0

    ϕ_x[:,2] = 0.5*ones(size(r))/c0
    ϕ_x[:,1] = -ϕ_x[:,2]

    djac = c0

    w_g = djac*w

  elseif (n==3)
    # Transform coordinates for quadratic elements
    c0 = ( x[n]-x[1] )/2.0
    c1 = ( x[n]+x[1] )/2.0
   
    x_g = c0.*r .+ c1

    ϕ[:,1] = 0.5.*r.*( r.-1.0 )
    ϕ[:,2] = -( r.+1.0 ).*( r.-1.0 )
    ϕ[:,3] = 0.5*r.*( r.+1.0 )

    ϕ_x[:,1] = ( r.-0.5 )/c0
    ϕ_x[:,2] = -2.0*r./c0
    ϕ_x[:,3] = ( r.+0.5 )./c0

    djac = c0

    w_g = djac*w

    ϕ_xx[:,1] = ones(size(r))/c0^2
    ϕ_xx[:,2] =-2.0*ϕ_xx[:,1]
    ϕ_xx[:,3] = ϕ_xx[:,1]

  elseif (n==4)
    # Transform coordinates for (nonconforming) cubic elements
    c0 = ( x[n]-x[1] )/2.0
    c1 = ( x[n]+x[1] )/2.0
    
    x_g = c0*r + c1
    
    r2  = r.*r
    r3  = r.*r2
 
    ϕ[:,1] = - 9.0*( r3-r2-r/9.0+1.0/9.0 )/16.0
    ϕ[:,2] =  27.0*( r3-r2/3.0-r+1.0/3.0 )/16.0
    ϕ[:,3] = -27.0*( r3+r2/3.0-r-1.0/3.0 )/16.0
    ϕ[:,4] =   9.0*( r3+r2-r/9.0-1.0/9.0 )/16.0
 
    ϕ_r[:,1] = - 9.0*( 3.0*r2-2.0*r-1.0/9.0 )/16.0
    ϕ_r[:,2] =  27.0*( 3.0*r2-2.0*r/3.0-1.0 )/16.0
    ϕ_r[:,3] = -27.0*( 3.0*r2+2.0*r/3.0-1.0 )/16.0
    ϕ_r[:,4] =   9.0*( 3.0*r2+2.0*r-1.0/9.0 )/16.0
 
    ϕ_rr[:,1] = - 9.0*( 6.0*r-2.0     )/16.0
    ϕ_rr[:,2] =  27.0*( 6.0*r-2.0/3.0 )/16.0
    ϕ_rr[:,3] = -27.0*( 6.0*r+2.0/3.0 )/16.0
    ϕ_rr[:,4] =   9.0*( 6.0*r+2.0     )/16.0
 
    dxdr = ϕ_r*x(:,1)
    djac = dxdr
    drdx = 1.0./djac
 
    ϕ_x[:,1] = ϕ_r[:,1].*drdx
    ϕ_x[:,2] = ϕ_r[:,2].*drdx
    ϕ_x[:,3] = ϕ_r[:,3].*drdx
    ϕ_x[:,4] = ϕ_r[:,4].*drdx
    w_g = djac.*w

  else
    error("Elements higher than cubic not currently supported")

  end

  return x_g, w_g, ϕ, ϕ_x, ϕ_xx

end
