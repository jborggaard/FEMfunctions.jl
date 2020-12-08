function onedQuadratureRule(rule::Int)
#= """
  onedQuadratureRule - Return Gauss integration points over (-1,1)

  Author: Jeff Borggaard, Virginia Tech
  Version: 1.3

  Usage:
  ```julia
    r,w = onedQuadratureRule(rule)
  ```

  Arguments:
  - `rule`: Number of Gauss points:

  Output arguments:
  - `r`: Gauss points located between (-1,1)      
  - `w`: Gauss weights corresponding to r
""" =#

  r = zeros(Float64,rule);
  w = zeros(Float64,rule);

  if rule == 1         # up to order 1 polynomials integrated exactly
    r[1] = 0.0
    w[1] = 2.0
    
  elseif rule == 2     # up to order 3 polynomials integrated exactly
    r[1] =-1.0 / sqrt(3.0);
    r[2] =-r[1]
    w[1] = 1.0
    w[2] = 1.0
    
  elseif rule == 3     # up to order 5 polynomials integrated exactly
    r[1] =-sqrt(3.0/5.0)
    r[2] = 0.0
    r[3] =-r[1]
    w[1] = 5.0 / 9.0
    w[2] = 8.0 / 9.0
    w[3] = w[1]
    
  elseif rule == 4    # up to order 7 polynomials integrated exactly
    r[1] =-sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0)
    r[2] =-sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0)
    r[3] =-r[2]
    r[4] =-r[1]
    w[1] = 0.5 - 1.0 / ( 6.0 * sqrt(6.0/5.0) )
    w[2] = 0.5 + 1.0 / ( 6.0 * sqrt(6.0/5.0) )
    w[3] = w[2]
    w[4] = w[1]
    
  elseif rule == 5    # up to order 9 polynomials integrated exactly
    r[1] =-sqrt(5.0+4.0*sqrt(5.0/14.0)) / 3.0
    r[2] =-sqrt(5.0-4.0*sqrt(5.0/14.0)) / 3.0
    r[3] = 0.0
    r[4] =-r[2]
    r[5] =-r[1]
    w[1] = 161.0/450.0-13.0/(180.0*sqrt(5.0/14.0))
    w[2] = 161.0/450.0+13.0/(180.0*sqrt(5.0/14.0))
    w[3] = 128.0/225.0
    w[4] = w[2]
    w[5] = w[1]
    
  elseif rule == 6    # up to order 11 polynomials integrated exactly
    r[1] = -0.2386191860831969
    r[2] = -0.6612093864662645
    r[3] = -0.9324695142031521
    r[4] = - r[1]
    r[5] = - r[2]
    r[6] = - r[3]
    w[1] = 0.4679139345726910
    w[2] = 0.3607615730481386
    w[3] = 0.1713244923791704
    w[4] = w[1]
    w[5] = w[2]
    w[6] = w[3]
    
  elseif (rule == 7)    # up to order 13 polynomials integrated exactly
    r[1] = -0.9491079123427585
    r[2] = -0.7415311855993945
    r[3] = -0.4058451513773972
    r[4] =  0.0000000000000000
    r[5] = - r[3]
    r[6] = - r[2]
    r[7] = - r[1]
    w[1] = 0.1294849661688697
    w[2] = 0.2797053914892766
    w[3] = 0.3818300505051189
    w[4] = 0.4179591836734694
    w[5] = w[3]
    w[6] = w[2]
    w[7] = w[1]
    
  elseif (rule == 8)    # up to order 15 polynomials integrated exactly
    r[1] = -0.9602898564975363
    r[2] = -0.7966664774136267
    r[3] = -0.5255324099163290
    r[4] = -0.1834346424956498
    r[5] = - r[4]
    r[6] = - r[3]
    r[7] = - r[2]
    r[8] = - r[1]
    w[1] = 0.1012285362903763
    w[2] = 0.2223810344533745
    w[3] = 0.3137066458778873
    w[4] = 0.3626837833783620
    w[5] = w[4]
    w[6] = w[3]
    w[7] = w[2]
    w[8] = w[1]

  elseif (rule == 9)     #  up to order 17 polynomials integrated exactly
    r[1] = -0.9681602395076261
    r[2] = -0.8360311073266358
    r[3] = -0.6133714327005904
    r[4] = -0.3242534234038089
    r[5] =  0.0000000000000000
    r[6] = - r[4]
    r[7] = - r[3]
    r[8] = - r[2]
    r[9] = - r[1]
    w[1] = 0.0812743883615744
    w[2] = 0.1806481606948574
    w[3] = 0.2606106964029354
    w[4] = 0.3123470770400029
    w[5] = 0.3302393550012598
    w[6] = w[4]
    w[7] = w[3]
    w[8] = w[2]
    w[9] = w[1]
  
  elseif (rule == 10)    #  up to order 19 polynomials integrated exactly
    r[ 1] = -0.9739065285171717
    r[ 2] = -0.8650633666889845
    r[ 3] = -0.6794095682990244
    r[ 4] = -0.4333953941292472
    r[ 5] = -0.1488743389816312
    r[ 6] = - r[5]
    r[ 7] = - r[4]
    r[ 8] = - r[3]
    r[ 9] = - r[2]
    r[10] = - r[1]
    w[ 1] = 0.0666713443086881
    w[ 2] = 0.1494513491505806
    w[ 3] = 0.2190863625159820
    w[ 4] = 0.2692667193099963
    w[ 5] = 0.2955242247147529
    w[ 6] = w[5]
    w[ 7] = w[4]
    w[ 8] = w[3]
    w[ 9] = w[2]
    w[10] = w[1]

  elseif (rule == 11)     #  up to order 21 integrated exactly
    r[ 1] = -0.9782286581460570
    r[ 2] = -0.8870625997680953
    r[ 3] = -0.7301520055740494
    r[ 4] = -0.5190961292068118
    r[ 5] = -0.2695431559523450
    r[ 6] =  0.0000000000000000
    r[ 7] = - r[5]
    r[ 8] = - r[4]
    r[ 9] = - r[3]
    r[10] = - r[2]
    r[11] = - r[1]
    w[ 1] = 0.0556685671161737
    w[ 2] = 0.1255803694649046
    w[ 3] = 0.1862902109277343
    w[ 4] = 0.2331937645919905
    w[ 5] = 0.2628045445102467
    w[ 6] = 0.2729250867779006
    w[ 7] = w[5]
    w[ 8] = w[4]
    w[ 9] = w[3]
    w[10] = w[2]
    w[11] = w[1]

  else
    error("onedQuadratureRule, rule not supported")
    
  end

  return r, w
end
