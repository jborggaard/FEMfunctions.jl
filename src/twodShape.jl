function twodShape(x,r,s,w)
#= """
  twodShape - computes test functions and derivatives on a triangular element
              given the element coordinates and quadrature points.

              ! Note: optimized for straight-sided elements.  Use
              ! `twod_shapeiso' for isoparametric elements.

              local unknown numbering follows the conventions

            3             3             3             3
            |\            |\            |\            87
            | \     or    6 5     or    6 5     or    | \   x is 10
            |  \          |  \          |7 \          9x 6
            1---2         1-4-2         1-4-2         14-52

             P1            P2      Crouzeux-Raviart    P3


  Author: Jeff Borggaard, Virginia Tech
  Version: 1.1

  Usage:
  ```julia
    xg,wg,ϕ,ϕ_x,ϕ_y = twodShape(x,r,s,w)
  ```

  Arguments:
  - `x`: Coordinates of the element nodes
  - `(r,s)`: Coordinates of Gauss points in unit triangle
  - `w`: Gauss weights associated with (r,s)

  Output arguments:
  - `xg`: Coordinates of Gauss points in the element
  - `wg`: Gauss weights scaled by the element Jacobian
  - `ϕ`: Value of element shape functions at xg
  - `ϕ_x`:
  - `ϕ_y`: First spatial derivatives of ϕ

  Changes:
      Crouzeux-Raviart element implemented by Alexander Hay, 2007
""" =#

  n   = size(x,1)
  ng  = size(r,1)

  xg  = zeros(Float64,ng,2)
  ϕ   = zeros(Float64,ng,n)
  ϕ_x = zeros(Float64,ng,n)
  ϕ_y = zeros(Float64,ng,n)

  wg  = zeros(Float64,ng)

  o   = ones(Float64,size(r))

  # Compute (r,s) -> (x,y) transformation for straight-sided elements
  c0 =  x[1,:]
  c1 = -x[1,:] + x[2,:]
  c2 = -x[1,:]          + x[3,:]

  xg[:,1] = c0[1]*o .+ c1[1]*r .+ c2[1]*s
  xr      = c1[1]*o
  xs      = c2[1]*o

  xg[:,2] = c0[2]*o .+ c1[2]*r .+ c2[2]*s
  yr      = c1[2]*o
  ys      = c2[2]*o

  # Compute the Jacobian of the (r,s) -> (x,y) transformation
  jac = xr.*ys .- yr.*xs
  wg  = jac.*w

  rx  = ys./jac
  sx  =-yr./jac
  ry  =-xs./jac
  sy  = xr./jac


  # Compute shape function and derivatives at Gauss points
  if n == 3
    ϕ[:,1] = 1.0*o - r  - s
    ϕ[:,2] =         r     
    ϕ[:,3] =              s
   
    ϕ_x[:,1] =       -rx - sx
    ϕ_x[:,2] =        rx     
    ϕ_x[:,3] =             sx
   
    ϕ_y[:,1] =      - ry - sy
    ϕ_y[:,2] =        ry     
    ϕ_y[:,3] =             sy

  elseif n == 6 
    ϕ[:,1] = 1.0*o - 3.0*r - 3.0*s + 2.0*r.*r + 4.0*r.*s + 2.0*s.*s
    ϕ[:,2] =       - 1.0*r         + 2.0*r.*r                      
    ϕ[:,3] =               - 1.0*s                       + 2.0*s.*s
    ϕ[:,4] =         4.0*r         - 4.0*r.*r - 4.0*r.*s           
    ϕ[:,5] =                                    4.0*r.*s           
    ϕ[:,6] =                 4.0*s            - 4.0*r.*s - 4.0*s.*s
  
    ϕ_x[:,1] = ( -3.0*o + 4.0*r + 4.0*s ).*rx + ( -3.0*o + 4.0*r + 4.0*s ).*sx
    ϕ_x[:,2] = ( -1.0*o + 4.0*r         ).*rx                                 
    ϕ_x[:,3] =                                  ( -1.0*o         + 4.0*s ).*sx
    ϕ_x[:,4] = (  4.0*o - 8.0*r - 4.0*s ).*rx + (        - 4.0*r         ).*sx
    ϕ_x[:,5] = (                  4.0*s ).*rx + (          4.0*r         ).*sx
    ϕ_x[:,6] = (                - 4.0*s ).*rx + (  4.0*o - 4.0*r - 8.0*s ).*sx
   
    ϕ_y[:,1] = ( -3.0*o + 4.0*r + 4.0*s ).*ry + ( -3.0*o + 4.0*r + 4.0*s ).*sy
    ϕ_y[:,2] = ( -1.0*o + 4.0*r         ).*ry                                 
    ϕ_y[:,3] =                                  ( -1.0*o         + 4.0*s ).*sy
    ϕ_y[:,4] = (  4.0*o - 8.0*r - 4.0*s ).*ry + (        - 4.0*r         ).*sy
    ϕ_y[:,5] = (                  4.0*s ).*ry + (          4.0*r         ).*sy
    ϕ_y[:,6] = (                - 4.0*s ).*ry + (  4.0*o - 4.0*r - 8.0*s ).*sy
    
  elseif n == 7
    ϕ[:,1] = (1.0*o-r-s).*(2.0*(1.0-r-s)-1.0) +  3.0*(1.0-r-s).*r.*s 
    ϕ[:,2] = r.*(2.0*r-1.0)                   +  3.0*(1.0-r-s).*r.*s 
    ϕ[:,3] = s.*(2.0*s-1.0)                   +  3.0*(1.0-r-s).*r.*s 
    ϕ[:,4] = 4.0*(1.0-r-s).*r                 - 12.0*(1.0-r-s).*r.*s 
    ϕ[:,5] = 4.0*r.*s                         - 12.0*(1.0-r-s).*r.*s 
    ϕ[:,6] = 4.0*s.*(1.0-r-s)                 - 12.0*(1.0-r-s).*r.*s 
    ϕ[:,7] = 27.0*(1.0-r-s).*r.*s                   

    ϕ_r = zeros(Float64,ng,n)
    ϕ_r[:,1] = -3.0 + 4.0*r + 7.0*s - 6.0*r.*s - 3.0*(s.^2)
    ϕ_r[:,2] = -1.0 + 4.0*r + 3.0*s - 6.0*r.*s - 3.0*(s.^2)
    ϕ_r[:,3] =                3.0*s - 6.0*r.*s - 3.0*(s.^2)
    ϕ_r[:,4] =  4.0 - 8.0*r -16.0*s +24.0*r.*s +12.0*(s.^2)
    ϕ_r[:,5] =              - 8.0*s +24.0*r.*s +12.0*(s.^2)
    ϕ_r[:,6] =              -16.0*s +24.0*r.*s +12.0*(s.^2)
    ϕ_r[:,7] =               27.0*s -54.0*r.*s -27.0*(s.^2)

    ϕ_s = zeros(Float64,ng,n)
    ϕ_s[:,1] = -3.0 + 7.0*r + 4.0*s - 6.0*r.*s - 3.0*(r.^2)
    ϕ_s[:,2] =        3.0*r         - 6.0*r.*s - 3.0*(r.^2)
    ϕ_s[:,3] = -1.0 + 3.0*r + 4.0*s - 6.0*r.*s - 3.0*(r.^2)
    ϕ_s[:,4] =      -16.0*r         +24.0*r.*s +12.0*(r.^2)
    ϕ_s[:,5] =      - 8.0*r         +24.0*r.*s +12.0*(r.^2)
    ϕ_s[:,6] =  4.0 -16.0*r - 8.0*s +24.0*r.*s +12.0*(r.^2)
    ϕ_s[:,7] =       27.0*r         -54.0*r.*s -27.0*(r.^2)

    ϕ_x = ϕ_r.*rx + ϕ_s.*sx
    ϕ_y = ϕ_r.*ry + ϕ_s.*sy

  elseif ( n==10 )
    # Compute shape function and derivatives at Gauss points
    t = ones(size(r)) - r - s
    
    ϕ = zeros(Float64,ng,n);
    ϕ[:, 1] =  4.5*(t-1.0/3.0).*(t-2.0/3.0).*t
    ϕ[:, 2] =  4.5*(r-1.0/3.0).*(r-2.0/3.0).*r
    ϕ[:, 3] =  4.5*(s-1.0/3.0).*(s-2.0/3.0).*s
    ϕ[:, 4] =  9.0*r-22.5*r.*(r+s)+13.5*r.*(r+s).^2
    ϕ[:, 5] = -4.5*r+18.0*r.^2+4.5*r.*s-13.5*r.^3-13.5*(r.^2).*s
    ϕ[:, 6] = -4.5*r.*s+13.5*(r.^2).*s
    ϕ[:, 7] = -4.5*r.*s+13.5*r.*s.^2
    ϕ[:, 8] = -4.5*s+4.5*r.*s+18.0*s.^2-13.5*r.*s.^2-13.5*s.^3
    ϕ[:, 9] =  9.0*s-22.5*r.*s-22.5*s.^2+13.5*s.*(r+s).^2
    ϕ[:,10] = 27.0*r.*s-27.0*(r.^2).*s-27.0*r.*s.^2

    ϕ_r = zeros(Float64,ng,n);
    ϕ_r[:, 1] = -5.5 +18.0*(r+s)-13.5*(r+s).^2
    ϕ_r[:, 2] =  1.0 -9.0*r+13.5*r.^2
#   ϕ_r[:, 3] =  zeros(Float64,ng,1)
    ϕ_r[:, 4] =  9.0 -45.0*r-22.5*s+40.5*r.^2+54.0*r.*s+13.5*s.^2
    ϕ_r[:, 5] = -4.5 +36.0*r+4.5*s-40.5*r.^2-27.0*r.*s
    ϕ_r[:, 6] = -4.5*s+27.0*r.*s
    ϕ_r[:, 7] = -4.5*s+13.5*s.*s
    ϕ_r[:, 8] =  4.5*s-13.5*s.*s
    ϕ_r[:, 9] =-22.5*s+27.0*(r.*s+s.*s)
    ϕ_r[:,10] = 27.0*(s - 2.0*r.*s-s.*s)

    ϕ_s = zeros(Float64,ng,n)
    ϕ_s[:, 1] = -5.5+18.0*(r+s)-13.5*(r+s).^2
#   ϕ_s[:, 2] = zeros(Float64,ng,1)
    ϕ_s[:, 3] = 1.0-9.0*s+13.5*s.^2
    ϕ_s[:, 4] = -22.5*r+27.0*(r.^2+r.*s)
    ϕ_s[:, 5] = 4.5*r-13.5*r.^2
    ϕ_s[:, 6] =-4.5*r+13.5*r.^2
    ϕ_s[:, 7] =-4.5*r+27.0*r.*s
    ϕ_s[:, 8] =-4.5+4.5*r+36.0*s-27.0*r.*s-40.5*s.^2
    ϕ_s[:, 9] = 9.0-22.5*r-45.0*s+13.5*r.^2+54.0*r.*s+40.5*s.^2
    ϕ_s[:,10] = 27.0*(r - r.^2 - 2.0*r.*s)

    ϕ_x = ϕ_r.*rx + ϕ_s.*sx
    ϕ_y = ϕ_r.*ry + ϕ_s.*sy
    
#   else
#     warning('element not supported')
  end

  return xg, wg, ϕ, ϕ_x, ϕ_y

end
