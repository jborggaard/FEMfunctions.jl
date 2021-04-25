function threedQuadratureRule(rule)
#=
   threed_quadrature.m - calculate quadrature integration points for 
                         tetrahedral elements
 
   Copyright (c) 2002, Jeff Borggaard, Virginia Tech
   Version: 1.0
 
   Usage:    [r,s,t,w] = threed_quadrature(rule)
 
   Variables:     rule
                         Number of quadrature points:
                            rule =  1, 1st degree monomial terms
                                    4, 2nd degree monomial terms
                                    5, 3rd degree monomial terms
                                   11, 4th degree monomial terms
                                   15, 5th degree monomial terms
                                   24, 6th degree monomial terms
                                   31, 7th degree monomial terms
                                   41, 8th degree monomial terms
 
                  r
                         xi coordinate of quadrature points
                  s
                         eta coordinate of quadrature points
                  t
                         zeta coordinate of quadrature points
                  w
                         quadrature weights corresponding to (r,s,t)
 
   Rules 11, 15, 24, 31, and 45 are from Keast, Moderate-degree 
   tetrahedral quadrature formulas, CMAME, v.55, (1986), 339-348.
%-----------------------------------------------------------------------
=#

  r = zeros(Float64,rule)
  s = zeros(Float64,rule)
  t = zeros(Float64,rule)
  w = zeros(Float64,rule)

  if (rule == 1)
    # The following points correspond to a 1 point rule
    # Monomials of degree 1 are integrated exactly
    r[1] = 0.250000000000000; s[1] = 0.250000000000000; t[1] = 0.250000000000000

    w[1] = 0.166666666666667  # (the volume of a regular tet is 1/6)

  elseif (rule == 4)
    # The following points correspond to a 4 point rule
    # Monomials of degree 2 are integrated exactly
    r[1] = 0.13819660;  s[1] = 0.13819660;  t[1] = 0.13819660
    r[2] = 0.58541020;  s[2] = 0.13819660;  t[2] = 0.13819660
    r[3] = 0.13819660;  s[3] = 0.58541020;  t[3] = 0.13819660
    r[4] = 0.13819660;  s[4] = 0.13819660;  t[4] = 0.58541020

    w[1] = 4.16666666666667e-2
    w[2] = 4.16666666666667e-2
    w[3] = 4.16666666666667e-2
    w[4] = 4.16666666666667e-2

  elseif (rule == 5)
    # The following points correspond to a 5 point rule
    # Monomials of degree 3 are integrated exactly
    r[1] = 0.250000000000000; s[1] = 0.250000000000000; t[1] = 0.250000000000000
    r[2] = 0.166666666666667; s[2] = 0.166666666666667; t[2] = 0.166666666666667
    r[3] = 0.500000000000000; s[3] = 0.166666666666667; t[3] = 0.166666666666667
    r[4] = 0.166666666666667; s[4] = 0.500000000000000; t[4] = 0.166666666666667
    r[5] = 0.166666666666667; s[5] = 0.166666666666667; t[5] = 0.500000000000000

    w[1] =-0.133333333333333
    w[2] = 0.075000000000000
    w[3] = 0.075000000000000
    w[4] = 0.075000000000000
    w[5] = 0.075000000000000

  elseif (rule == 11)
    # The following points correspond to an 11 point rule
    # Monomials of degree 4 are integrated exactly
    r[1]  = 0.25;          s[1]  = 0.25;          t[1]  = 0.25

    t1 = 0.0714285714285714285714285
    t2 = 0.785714285714285714285714

    r[2]  = t1;            s[2]  = t1;            t[2]  = t1
    r[3]  = t2;            s[3]  = t1;            t[3]  = t1
    r[4]  = t1;            s[4]  = t2;            t[4]  = t1
    r[5]  = t1;            s[5]  = t1;            t[5]  = t2

    t1 = 0.399403576166799219
    t2 = 0.100596423833200785

    r[6]  = t1;            s[6]  = t1;            t[6]  = t2
    r[7]  = t2;            s[7]  = t1;            t[7]  = t1
    r[8]  = t2;            s[8]  = t2;            t[8]  = t1
    r[9]  = t1;            s[9]  = t2;            t[9]  = t2
    r[10] = t2;            s[10] = t1;            t[10] = t2
    r[11] = t1;            s[11] = t2;            t[11] = t1

    w[1]  =-1.31555555555555550e-2

    w[2]  = 7.62222222222222222e-3
    w[3]  = 7.62222222222222222e-3
    w[4]  = 7.62222222222222222e-3
    w[5]  = 7.62222222222222222e-3

    w[6]  = 2.48888888888888880e-2
    w[7]  = 2.48888888888888880e-2
    w[8]  = 2.48888888888888880e-2
    w[9]  = 2.48888888888888880e-2
    w[10] = 2.48888888888888880e-2
    w[11] = 2.48888888888888880e-2
    
  elseif (rule == 15)
    # The following points correspond to a 15 point rule
    # Monomials of degree 5 are integrated exactly
    r[1]  = 0.25;          s[1]  = 0.25;          t[1]  = 0.25

    t1 = 0.333333333333333333333333
    t2 = 0.000000000000000000000000

    r[2]  = t1;            s[2]  = t1;            t[2]  = t1  
    r[3]  = t2;            s[3]  = t1;            t[3]  = t1
    r[4]  = t1;            s[4]  = t2;            t[4]  = t1
    r[5]  = t1;            s[5]  = t1;            t[5]  = t2

    t1 = 0.0909090909090909090909
    t2 = 0.7272727272727272727272

    r[6]  = t1;            s[6]  = t1;            t[6]  = t1
    r[7]  = t2;            s[7]  = t1;            t[7]  = t1
    r[8]  = t1;            s[8]  = t2;            t[8]  = t1
    r[9]  = t1;            s[9]  = t1;            t[9]  = t2

    t1 = 0.0665501535736642813
    t2 = 0.433449846426335728

    r[10] = t2;            s[10] = t2;            t[10] = t1
    r[11] = t1;            s[11] = t2;            t[11] = t2
    r[12] = t1;            s[12] = t1;            t[12] = t2
    r[13] = t2;            s[13] = t1;            t[13] = t1
    r[14] = t1;            s[14] = t2;            t[14] = t1
    r[15] = t2;            s[15] = t1;            t[15] = t2

    w[1]  = 0.0302836780970891856

    w[2]  = 0.00602678571428571597
    w[3]  = 0.00602678571428571597
    w[4]  = 0.00602678571428571597
    w[5]  = 0.00602678571428571597

    w[6]  = 0.0116452490860289742
    w[7]  = 0.0116452490860289742
    w[8]  = 0.0116452490860289742
    w[9]  = 0.0116452490860289742

    w[10] = 0.0109491415613864534
    w[11] = 0.0109491415613864534
    w[12] = 0.0109491415613864534
    w[13] = 0.0109491415613864534
    w[14] = 0.0109491415613864534
    w[15] = 0.0109491415613864534

  elseif (rule == 24)
    # The following points correspond to a 24 point rule
    # Monomials of degree 6 are integrated exactly
    
    t1 = 0.214602871259151684
    t2 = 0.356191386222544953

    r[1]  = t1;      s[1]  = t1;      t[1]  = t1

    r[2]  = t2;      s[2]  = t1;      t[2]  = t1
    r[3]  = t1;      s[3]  = t2;      t[3]  = t1
    r[4]  = t1;      s[4]  = t1;      t[4]  = t2

    t1 = 0.0406739585346113397
    t2 = 0.877978124396165982

    r[5]  = t1;      s[5]  = t1;      t[5]  = t1

    r[6]  = t2;      s[6]  = t1;      t[6]  = t1
    r[7]  = t1;      s[7]  = t2;      t[7]  = t1
    r[8]  = t1;      s[8]  = t1;      t[8]  = t2

    t1 = 0.322337890142275646
    t2 = 0.0329863295731730594

    r[9]  = t1;      s[9]  = t1;      t[9]  = t1

    r[10] = t2;      s[10] = t1;      t[10] = t1
    r[11] = t1;      s[11] = t2;      t[11] = t1
    r[12] = t1;      s[12] = t1;      t[12] = t2

    t1 = 0.0636610018750175299
    t2 = 0.269672331458315867
    t3 = 0.603005664791649076

    r[13] = t1;      s[13] = t1;      t[13] = t2
    r[14] = t1;      s[14] = t1;      t[14] = t3
    r[15] = t1;      s[15] = t2;      t[15] = t1
    r[16] = t1;      s[16] = t2;      t[16] = t3
    r[17] = t1;      s[17] = t3;      t[17] = t1
    r[18] = t1;      s[18] = t3;      t[18] = t2
    r[19] = t2;      s[19] = t1;      t[19] = t1
    r[20] = t2;      s[20] = t1;      t[20] = t3
    r[21] = t2;      s[21] = t3;      t[21] = t1
    r[22] = t3;      s[22] = t1;      t[22] = t1
    r[23] = t3;      s[23] = t1;      t[23] = t2
    r[24] = t3;      s[24] = t2;      t[24] = t1

    w[1]  = 6.65379170969464506e-3
    w[2]  = 6.65379170969464506e-3
    w[3]  = 6.65379170969464506e-3
    w[4]  = 6.65379170969464506e-3

    w[5]  = 1.67953517588677620e-3
    w[6]  = 1.67953517588677620e-3
    w[7]  = 1.67953517588677620e-3
    w[8]  = 1.67953517588677620e-3

    w[9]  = 9.22619692394239843e-3
    w[10] = 9.22619692394239843e-3
    w[11] = 9.22619692394239843e-3
    w[12] = 9.22619692394239843e-3

    w[13] = 8.03571428571428248e-3
    w[14] = 8.03571428571428248e-3
    w[15] = 8.03571428571428248e-3
    w[16] = 8.03571428571428248e-3
    w[17] = 8.03571428571428248e-3
    w[18] = 8.03571428571428248e-3
    w[19] = 8.03571428571428248e-3
    w[20] = 8.03571428571428248e-3
    w[21] = 8.03571428571428248e-3
    w[22] = 8.03571428571428248e-3
    w[23] = 8.03571428571428248e-3
    w[24] = 8.03571428571428248e-3

  elseif (rule == 31)
    # The following points correspond to a 31 point rule
    # Monomials of degree 7 are integrated exactly
    
    t1 = 0.0782131923303186549
    t2 = 0.765360423009044044

    r[1]  = t1;      s[1]  = t1;      t[1]  = t1
    r[2]  = t2;      s[2]  = t1;      t[2]  = t1
    r[3]  = t1;      s[3]  = t2;      t[3]  = t1
    r[4]  = t1;      s[4]  = t1;      t[4]  = t2

    t1 = 0.121843216663904411
    t2 = 0.634470350008286765

    r[5]  = t1;      s[5]  = t1;      t[5]  = t1
    r[6]  = t2;      s[6]  = t1;      t[6]  = t1
    r[7]  = t1;      s[7]  = t2;      t[7]  = t1
    r[8]  = t1;      s[8]  = t1;      t[8]  = t2

    t1 = 0.332539164446420554
    t2 = 0.00238250666073834549

    r[9]  = t1;      s[9]  = t1;      t[9]  = t1
    r[10] = t2;      s[10] = t1;      t[10] = t1
    r[11] = t1;      s[11] = t2;      t[11] = t1
    r[12] = t1;      s[12] = t1;      t[12] = t2

    t1 = 0.1
    t2 = 0.2
    t3 = 0.6

    r[13] = t1;      s[13] = t1;      t[13] = t2
    r[14] = t1;      s[14] = t1;      t[14] = t3
    r[15] = t1;      s[15] = t2;      t[15] = t1
    r[16] = t1;      s[16] = t2;      t[16] = t3
    r[17] = t1;      s[17] = t3;      t[17] = t1
    r[18] = t1;      s[18] = t3;      t[18] = t2
    r[19] = t2;      s[19] = t1;      t[19] = t1
    r[20] = t2;      s[20] = t1;      t[20] = t3
    r[21] = t2;      s[21] = t3;      t[21] = t1
    r[22] = t3;      s[22] = t1;      t[22] = t1
    r[23] = t3;      s[23] = t1;      t[23] = t2
    r[24] = t3;      s[24] = t2;      t[24] = t1

    t1 = 0.0
    t2 = 0.5

    r[25] = t2;      s[25] = t2;      t[25] = t1
    r[26] = t1;      s[26] = t2;      t[26] = t2
    r[27] = t1;      s[27] = t1;      t[27] = t2
    r[28] = t2;      s[28] = t1;      t[28] = t1
    r[29] = t1;      s[29] = t2;      t[29] = t1
    r[30] = t2;      s[30] = t1;      t[30] = t2

    t1 = 0.25

    r[31] = t1;      s[31] = t1;      t[31] = t1

    w[1]  = 0.01059994152441416093
    w[2]  = 0.01059994152441416093
    w[3]  = 0.01059994152441416093
    w[4]  = 0.01059994152441416093

    w[5]  =-0.0625177401143299494
    w[6]  =-0.0625177401143299494
    w[7]  =-0.0625177401143299494
    w[8]  =-0.0625177401143299494

    w[9]  = 0.00489142526307353653
    w[10] = 0.00489142526307353653
    w[11] = 0.00489142526307353653
    w[12] = 0.00489142526307353653

    w[13] = 0.0275573192239850917
    w[14] = 0.0275573192239850917
    w[15] = 0.0275573192239850917
    w[16] = 0.0275573192239850917
    w[17] = 0.0275573192239850917
    w[18] = 0.0275573192239850917
    w[19] = 0.0275573192239850917
    w[20] = 0.0275573192239850917
    w[21] = 0.0275573192239850917
    w[22] = 0.0275573192239850917
    w[23] = 0.0275573192239850917
    w[24] = 0.0275573192239850917

    w[25] = 0.000970017636684296702
    w[26] = 0.000970017636684296702
    w[27] = 0.000970017636684296702
    w[28] = 0.000970017636684296702
    w[29] = 0.000970017636684296702
    w[30] = 0.000970017636684296702

    w[31] = 0.0182642234661087939

  elseif (rule == 45)
    # The following points correspond to a 45 point rule
    # Monomials of degree 8 are integrated exactly
    
    r[1]  = 0.25;    s[1]  = 0.25;    t[1]  = 0.25
    
    t1 = 0.127470936566639015
    t2 = 0.617587190300082967

    r[2]  = t1;      s[2]  = t1;      t[2]  = t1
    r[3]  = t2;      s[3]  = t1;      t[3]  = t1
    r[4]  = t1;      s[4]  = t2;      t[4]  = t1
    r[5]  = t1;      s[5]  = t1;      t[5]  = t2

    t1 = 0.0320788303926322960
    t2 = 0.903763508822103123

    r[6]  = t1;      s[6]  = t1;      t[6]  = t1
    r[7]  = t2;      s[7]  = t1;      t[7]  = t1
    r[8]  = t1;      s[8]  = t2;      t[8]  = t1
    r[9]  = t1;      s[9]  = t1;      t[9]  = t2

    t1 = 0.0497770956432810185
    t2 = 0.450222904356718978

    r[10] = t2;      s[10] = t2;      t[10] = t1
    r[11] = t1;      s[11] = t2;      t[11] = t2
    r[12] = t1;      s[12] = t1;      t[12] = t2
    r[13] = t2;      s[13] = t1;      t[13] = t1
    r[14] = t1;      s[14] = t2;      t[14] = t1
    r[15] = t2;      s[15] = t1;      t[15] = t2

    t1 = 0.183730447398549945
    t2 = 0.316269552601450060

    r[16] = t2;      s[16] = t2;      t[16] = t1
    r[17] = t1;      s[17] = t2;      t[17] = t2
    r[18] = t1;      s[18] = t1;      t[18] = t2
    r[19] = t2;      s[19] = t1;      t[19] = t1
    r[20] = t1;      s[20] = t2;      t[20] = t1
    r[21] = t2;      s[21] = t1;      t[21] = t2

    t1 = 0.231901089397150906
    t2 = 0.0229177878448171174
    t3 = 0.513280033360881072

    r[22] = t1;      s[22] = t1;      t[22] = t2
    r[23] = t1;      s[23] = t1;      t[23] = t3
    r[24] = t1;      s[24] = t2;      t[24] = t1
    r[25] = t1;      s[25] = t2;      t[25] = t3
    r[26] = t1;      s[26] = t3;      t[26] = t1
    r[27] = t1;      s[27] = t3;      t[27] = t2
    r[28] = t2;      s[28] = t1;      t[28] = t1
    r[29] = t2;      s[29] = t1;      t[29] = t3
    r[30] = t2;      s[30] = t3;      t[30] = t1
    r[31] = t3;      s[31] = t1;      t[31] = t1
    r[32] = t3;      s[32] = t1;      t[32] = t2
    r[33] = t3;      s[33] = t2;      t[33] = t1

    t1 = 0.0379700484718286102
    t2 = 0.730313427807538396
    t3 = 0.193746475248804382

    r[34] = t1;      s[34] = t1;      t[34] = t2
    r[35] = t1;      s[35] = t1;      t[35] = t3
    r[36] = t1;      s[36] = t2;      t[36] = t1
    r[37] = t1;      s[37] = t2;      t[37] = t3
    r[38] = t1;      s[38] = t3;      t[38] = t1
    r[39] = t1;      s[39] = t3;      t[39] = t2
    r[40] = t2;      s[40] = t1;      t[40] = t1
    r[41] = t2;      s[41] = t1;      t[41] = t3
    r[42] = t2;      s[42] = t3;      t[42] = t1
    r[43] = t3;      s[43] = t1;      t[43] = t1
    r[44] = t3;      s[44] = t1;      t[44] = t2
    r[45] = t3;      s[45] = t2;      t[45] = t1

    w[1]  =-0.0393270066412926145

    w[2]  = 0.00408131605934270525
    w[3]  = 0.00408131605934270525
    w[4]  = 0.00408131605934270525
    w[5]  = 0.00408131605934270525

    w[6]  = 0.000658086773304341943
    w[7]  = 0.000658086773304341943
    w[8]  = 0.000658086773304341943
    w[9]  = 0.000658086773304341943

    w[10] = 0.00438425882512284693
    w[11] = 0.00438425882512284693
    w[12] = 0.00438425882512284693
    w[13] = 0.00438425882512284693
    w[14] = 0.00438425882512284693
    w[15] = 0.00438425882512284693

    w[16] = 0.0138300638425098166
    w[17] = 0.0138300638425098166
    w[18] = 0.0138300638425098166
    w[19] = 0.0138300638425098166
    w[20] = 0.0138300638425098166
    w[21] = 0.0138300638425098166

    w[22] = 0.00424043742468372453
    w[23] = 0.00424043742468372453
    w[24] = 0.00424043742468372453
    w[25] = 0.00424043742468372453
    w[26] = 0.00424043742468372453
    w[27] = 0.00424043742468372453
    w[28] = 0.00424043742468372453
    w[29] = 0.00424043742468372453
    w[30] = 0.00424043742468372453
    w[31] = 0.00424043742468372453
    w[32] = 0.00424043742468372453
    w[33] = 0.00424043742468372453

    w[34] = 0.00223873973961420164
    w[35] = 0.00223873973961420164
    w[36] = 0.00223873973961420164
    w[37] = 0.00223873973961420164
    w[38] = 0.00223873973961420164
    w[39] = 0.00223873973961420164
    w[40] = 0.00223873973961420164
    w[41] = 0.00223873973961420164
    w[42] = 0.00223873973961420164
    w[43] = 0.00223873973961420164
    w[44] = 0.00223873973961420164
    w[45] = 0.00223873973961420164

  else
      r=0
      s=0
      t=0
      w=0
    #error('quadrature rules other than 1, 4, 5, 11, 15, 24, 31, or 45 are not supported\n')
  end

  return r, s, t, w
end
