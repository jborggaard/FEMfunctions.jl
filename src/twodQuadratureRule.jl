function twodQuadratureRule(rule::Int)
#= """
  twodQuadratureRule - calculate numerical integration points for triangular
                       elements

  Author: Jeff Borggaard, Virginia Tech
          part of FEMfunctions

  Licensing:
     This code is distributed under the MIT license.
 
  Usage:
  ```julia
    r,s,w = twodQuadratureRule(rule)
  ```

  Arguments:
  - `r`: xi coordinate of Gauss points
  - `s`: eta coordinate of Gauss points
  - `w`: Gauss weights corresponding to (r,s)

  Output arguments:
  - `rule`: Number of Gauss points:
""" =#

  r = zeros(Float64,rule)
  s = zeros(Float64,rule)
  w = zeros(Float64,rule)

  if rule == 1
    # The trivial linear triangle case
    r[1] = 1.0/3.0       
    s[1] = 1.0/3.0

    w[1] = 0.5;

  elseif rule == 3
    # The following points correspond to a 3 point rule (quadratics)
    r[1] = 2.0/3.0;         s[1] = 1.0/6.0
    r[2] = 1.0/6.0;         s[2] = 2.0/3.0
    r[3] = 1.0/6.0;         s[3] = 1.0/6.0

    w[1] = 1.0/6.0
    w[2] = w[1]
    w[3] = w[1]

  elseif rule == 7
    # The following points correspond to a 7 point rule,
    # see Dunavant, IJNME, v. 21, pp. 1129-1148, 1995.
    # or Braess, p. 95.   (integrates 5th degree)

    t1 = 1.0/3.0;  t2 = (6.0+sqrt(15.0))/21.0;   t3 = 4.0/7.0 - t2

    r[1] = t1;            s[1] = t1
    r[2] = t2;            s[2] = t2
    r[3] = 1.0-2.0*t2;    s[3] = t2
    r[4] = t2;            s[4] = r[3]
    r[5] = t3;            s[5] = t3
    r[6] = 1.0-2.0*t3;    s[6] = t3
    r[7] = t3;            s[7] = r[6]

    t1 = 9.0/80.0;  t2 = (155.0+sqrt(15.0))/2400.0;  t3 = 31.0/240.0 - t2

    w[1]  = t1
    w[2]  = t2
    w[3]  = t2
    w[4]  = t2
    w[5]  = t3
    w[6]  = t3
    w[7]  = t3

  elseif rule == 13   # same Dunavant ref.  (integrates 7th degree)
    r[1]  = 0.0651301029022;  s[1]  = r[1]
    r[2]  = 0.8697397941956;  s[2]  = r[1]
    r[3]  = r[1];             s[3]  = r[2]
    r[4]  = 0.3128654960049;  s[4]  = 0.0486903154253
    r[5]  = 0.6384441885698;  s[5]  = r[4]
    r[6]  = s[4];             s[6]  = r[5]
    r[7]  = r[5];             s[7]  = r[6]
    r[8]  = r[4];             s[8]  = r[5]
    r[9]  = r[6];             s[9]  = r[4]
    r[10] = 0.2603459660790;  s[10] = r[10]
    r[11] = 0.4793080678419;  s[11] = r[10]
    r[12] = r[10];            s[12] = r[11]
    r[13] = 0.3333333333333;  s[13] = r[13]

    w[1]  = 0.0533472356088
    w[2]  = w[1]
    w[3]  = w[1]
    w[4]  = 0.0771137608903
    w[5]  = w[4]
    w[6]  = w[4]
    w[7]  = w[4]
    w[8]  = w[4]
    w[9]  = w[4]
    w[10] = 0.1756152574332
    w[11] = w[10]
    w[12] = w[10]
    w[13] =-0.1495700444677

    w = w/2.0

  elseif rule == 19   # same Dunavant ref.  (integrates 9th degree)
    r[1]  = 0.036838412054736; s[1]  = 0.741198598784498
    r[2]  = 0.741198598784498; s[2]  = 0.221962989160766
    r[3]  = 0.221962989160766; s[3]  = 0.036838412054736
    r[4]  = 0.741198598784498; s[4]  = 0.036838412054736
    r[5]  = 0.221962989160766; s[5]  = 0.741198598784498
    r[6]  = 0.036838412054736; s[6]  = 0.221962989160766
    r[7]  = 0.044729513394453; s[7]  = 0.910540973211095
    r[8]  = 0.044729513394453; s[8]  = 0.044729513394453
    r[9]  = 0.910540973211095; s[9]  = 0.044729513394453
    r[10] = 0.188203535619033; s[10] = 0.623592928761935
    r[11] = 0.188203535619033; s[11] = 0.188203535619033
    r[12] = 0.623592928761935; s[12] = 0.188203535619033
    r[13] = 0.437089591492937; s[13] = 0.125820817014127
    r[14] = 0.437089591492937; s[14] = 0.437089591492937
    r[15] = 0.125820817014127; s[15] = 0.437089591492937
    r[16] = 0.489682519198738; s[16] = 0.020634961602525
    r[17] = 0.489682519198738; s[17] = 0.489682519198738
    r[18] = 0.020634961602525; s[18] = 0.489682519198738
    r[19] = 0.333333333333333; s[19] = 0.333333333333333

    w[1]  = 0.043283539377289
    w[2]  = w[1]
    w[3]  = w[1]
    w[4]  = w[1]
    w[5]  = w[1]
    w[6]  = w[1]
    w[7]  = 0.025577675658698
    w[8]  = w[7]
    w[9]  = w[7]
    w[10] = 0.079647738927210
    w[11] = w[10]
    w[12] = w[10]
    w[13] = 0.077827541004774
    w[14] = w[13]
    w[15] = w[14]
    w[16] = 0.031334700227139
    w[17] = w[16]
    w[18] = w[16]
    w[19] = 0.097135796282799

    w = w/2.0

  else
    error("No quadrature rules other than 1, 3, 7, 13 or 19 implemented")
  end

  return r, s, w

end
