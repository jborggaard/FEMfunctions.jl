function TetMesh_CubeMesh(xMin::Float64, xMax::Float64,
                          yMin::Float64, yMax::Float64,
                          zMin::Float64, zMax::Float64,
                          etype::String,
                          nNodesX::Int64, nNodesY::Int64, nNodesZ::Int64)
#= """
TetMesh_CubeMesh - Generate a mesh for a regular hexahedral domain (cube).  
                   Either linear or quadratic tetrahedral elements are
                   constructed from mini-"cubes."
 
   Author: Jeff Borggaard, Virginia Tech
           part of FEMfunctions.jl

   Licensing:
      This code is distributed under the MIT license.
 
   Usage:
   ```julia
     x, eConn = TetMesh_CubeMesh(xMin,xMax,yMin,yMax,zMin,zMax,
                                 etype,
                                 nNodesX,nNodesY,nNodesZ)
   ```

   Arguments:
   - `(xMin,xMax)`: range in x-direction
   - `(yMin,yMax)`: range in y-direction
   - `(zMin,zMax)`: range in z-direction

   - etype: Element type, 'linear' or 'quadratic'

   - `nNodesX`: number of nodes in the x direction
   - `nNodesY`: number of nodes in the y direction
   - `nNodesZ`: number of nodes in the z direction
  
   The number of nodes must be compatible with the element type,
   for quadratic elements, this means nNodes* is 4k+1 for natural number k

   Output arguments:
   - `x`: Nodal coordinates of the mesh
   - `eConn`: Element connectivity
""" =#

  #  Generate node coordinates
  nNodes = nNodesX*nNodesY*nNodesZ
  x = zeros(Float64,nNodes,3)
  dx = (xMax-xMin)/(nNodesX-1)  
  dy = (yMax-yMin)/(nNodesY-1)
  dz = (zMax-zMin)/(nNodesZ-1)
  for j=1:nNodesX
    pointsX[j] = dx*(j-1)
  end

  for k=1:nNodesY
    pointsY[k] = dy*(k-1)
  end

  for l=1:nNodesZ
    pointsZ[l] = dz*(l-1)
  end


  if etype == "linear"
    # preallocate storage
    x = zeros(nNodes,3)
    
    for l=1:nNodesZ
      for k=1:nNodesY
        for j=1:nNodesX
          i = j + (k-1)*nNodesX + (l-1)*nNodesX*nNodesY
          x[i,1] = pointsX[j]
          x[i,2] = pointsY[k]
          x[i,3] = pointsZ[l]
        end
      end
    end
  
  elseif etype == "quadratic"
    # preallocate storage
    nNodes = jkl_to_global(nNodesX-1,nNodesY-1,nNodesZ-1,nNodesX,nNodesY,lmax)
    x = zeros(nNodes,3)
    
    # borrow logic from C code
    for l=0:nNodesZ-1
      l_odd = mod(l,2)

      for k=0:nNodesY-1
        k_odd = mod(k,2)

        if k_odd && l_odd       # skip center of cubes
          for j=0:2:nNodesX-1
            i = jkl_to_global(j,k,l,nNodesX,nNodesY,nNodesZ)
            x[i,1] = pointsX[j+1]
            x[i,2] = pointsY[k+1]
            x[i,3] = pointsZ[l+1]
          end
        else                    # fill the entire row
          for j=0:nNodesX-1
            i = jkl_to_global(j,k,l,nNodesX,nNodesY,nNodesZ)
            x[i,1] = pointsX[j+1]
            x[i,2] = pointsY[k+1]
            x[i,3] = pointsZ[l+1]
          end
        end
      end # k-loop
    end # l-loop

  else
    error("TetMesh_cubeMesh: Orders other than 1 or 2 are not implemented.")
  end

  
  if etype=="linear"
    # Set up connectivity (for a mesh of linear elements)
    
    # preallocate storage
    nElems = (nNodesX-1)*(nNodesY-1)*(nNodesZ-1)*5
    eConn  = zeros(nElems,4)
    
    nel = 0;
    for l=1:nNodesZ-1
      for k=1:nNodesY-1
        for j=1:nNodesX-1
          i1 = j + (k-1)*nNodesX + (l-1)*nNodesX*nNodesY
          i2 = i1 + 1
          i3 = i1 + nNodesX
          i4 = i3 + 1
          i5 = i1 + nNodesX*nNodesY
          i6 = i5 + 1
          i7 = i5 + nNodesX
          i8 = i7 + 1

          eConn[nel+1,:] = [ i1 i5 i2 i3 ]
          eConn[nel+2,:] = [ i6 i2 i5 i8 ]
          eConn[nel+3,:] = [ i7 i3 i8 i5 ]
          eConn[nel+4,:] = [ i4 i8 i3 i2 ]
          eConn[nel+5,:] = [ i5 i8 i2 i3 ]  # connect the stump element

          warning("wrong_edges?, need to test")
          nel = nel + 5
        end
      end
    end

  elseif etype=="quadratic"
    # Set up connectivity (for a mesh of quadratic elements)
    # loop over 5x5x5 blocks and carve into 40 elements so that they tile...

    vertex_corner = false
    if vertex_corner  #  build a mesh with one vertex in the corners
      ee = 0
      for l=0:4:nNodesZ-3
      for k=0:4:nNodesY-3
      for j=0:4:nNodesX-3
        i0   = jkl_to_global(j  , k  , l  ,nNodesX,nNodesY,nNodesZ)
        i1   = jkl_to_global(j+1, k  , l  ,nNodesX,nNodesY,nNodesZ)
        i2   = jkl_to_global(j+2, k  , l  ,nNodesX,nNodesY,nNodesZ)
        i3   = jkl_to_global(j+3, k  , l  ,nNodesX,nNodesY,nNodesZ)
        i4   = jkl_to_global(j+4, k  , l  ,nNodesX,nNodesY,nNodesZ)
        i5   = jkl_to_global(j  , k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i6   = jkl_to_global(j+1, k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i7   = jkl_to_global(j+2, k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i8   = jkl_to_global(j+3, k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i9   = jkl_to_global(j+4, k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i10  = jkl_to_global(j  , k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i11  = jkl_to_global(j+1, k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i12  = jkl_to_global(j+2, k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i13  = jkl_to_global(j+3, k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i14  = jkl_to_global(j+4, k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i15  = jkl_to_global(j  , k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i16  = jkl_to_global(j+1, k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i17  = jkl_to_global(j+2, k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i18  = jkl_to_global(j+3, k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i19  = jkl_to_global(j+4, k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i20  = jkl_to_global(j  , k+4, l  ,nNodesX,nNodesY,nNodesZ)
        i21  = jkl_to_global(j+1, k+4, l  ,nNodesX,nNodesY,nNodesZ)
        i22  = jkl_to_global(j+2, k+4, l  ,nNodesX,nNodesY,nNodesZ)
        i23  = jkl_to_global(j+3, k+4, l  ,nNodesX,nNodesY,nNodesZ)
        i24  = jkl_to_global(j+4, k+4, l  ,nNodesX,nNodesY,nNodesZ)

        i25  = jkl_to_global(j  , k  , l+1,nNodesX,nNodesY,nNodesZ)
        i26  = jkl_to_global(j+1, k  , l+1,nNodesX,nNodesY,nNodesZ)
        i27  = jkl_to_global(j+2, k  , l+1,nNodesX,nNodesY,nNodesZ)
        i28  = jkl_to_global(j+3, k  , l+1,nNodesX,nNodesY,nNodesZ)
        i29  = jkl_to_global(j+4, k  , l+1,nNodesX,nNodesY,nNodesZ)
        i30  = jkl_to_global(j  , k+1, l+1,nNodesX,nNodesY,nNodesZ)
        i31  = jkl_to_global(j+1, k+1, l+1,nNodesX,nNodesY,nNodesZ)
        i32  = jkl_to_global(j+2, k+1, l+1,nNodesX,nNodesY,nNodesZ)
        i33  = jkl_to_global(j+3, k+1, l+1,nNodesX,nNodesY,nNodesZ)
        i34  = jkl_to_global(j+4, k+1, l+1,nNodesX,nNodesY,nNodesZ)
        i35  = jkl_to_global(j  , k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i36  = jkl_to_global(j+1, k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i37  = jkl_to_global(j+2, k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i38  = jkl_to_global(j+3, k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i39  = jkl_to_global(j+4, k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i40  = jkl_to_global(j  , k+3, l+1,nNodesX,nNodesY,nNodesZ)
        i41  = jkl_to_global(j+1, k+3, l+1,nNodesX,nNodesY,nNodesZ)
        i42  = jkl_to_global(j+2, k+3, l+1,nNodesX,nNodesY,nNodesZ)
        i43  = jkl_to_global(j+3, k+3, l+1,nNodesX,nNodesY,nNodesZ)
        i44  = jkl_to_global(j+4, k+3, l+1,nNodesX,nNodesY,nNodesZ)
        i45  = jkl_to_global(j  , k+4, l+1,nNodesX,nNodesY,nNodesZ)
        i46  = jkl_to_global(j+1, k+4, l+1,nNodesX,nNodesY,nNodesZ)
        i47  = jkl_to_global(j+2, k+4, l+1,nNodesX,nNodesY,nNodesZ)
        i48  = jkl_to_global(j+3, k+4, l+1,nNodesX,nNodesY,nNodesZ)
        i49  = jkl_to_global(j+4, k+4, l+1,nNodesX,nNodesY,nNodesZ)

        i50  = jkl_to_global(j  , k  , l+2,nNodesX,nNodesY,nNodesZ)
        i51  = jkl_to_global(j+1, k  , l+2,nNodesX,nNodesY,nNodesZ)
        i52  = jkl_to_global(j+2, k  , l+2,nNodesX,nNodesY,nNodesZ)
        i53  = jkl_to_global(j+3, k  , l+2,nNodesX,nNodesY,nNodesZ)
        i54  = jkl_to_global(j+4, k  , l+2,nNodesX,nNodesY,nNodesZ)
        i55  = jkl_to_global(j  , k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i56  = jkl_to_global(j+1, k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i57  = jkl_to_global(j+2, k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i58  = jkl_to_global(j+3, k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i59  = jkl_to_global(j+4, k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i60  = jkl_to_global(j  , k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i61  = jkl_to_global(j+1, k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i62  = jkl_to_global(j+2, k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i63  = jkl_to_global(j+3, k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i64  = jkl_to_global(j+4, k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i65  = jkl_to_global(j  , k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i66  = jkl_to_global(j+1, k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i67  = jkl_to_global(j+2, k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i68  = jkl_to_global(j+3, k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i69  = jkl_to_global(j+4, k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i70  = jkl_to_global(j  , k+4, l+2,nNodesX,nNodesY,nNodesZ)
        i71  = jkl_to_global(j+1, k+4, l+2,nNodesX,nNodesY,nNodesZ)
        i72  = jkl_to_global(j+2, k+4, l+2,nNodesX,nNodesY,nNodesZ)
        i73  = jkl_to_global(j+3, k+4, l+2,nNodesX,nNodesY,nNodesZ)
        i74  = jkl_to_global(j+4, k+4, l+2,nNodesX,nNodesY,nNodesZ)

        i75  = jkl_to_global(j  , k  , l+3,nNodesX,nNodesY,nNodesZ)
        i76  = jkl_to_global(j+1, k  , l+3,nNodesX,nNodesY,nNodesZ)
        i77  = jkl_to_global(j+2, k  , l+3,nNodesX,nNodesY,nNodesZ)
        i78  = jkl_to_global(j+3, k  , l+3,nNodesX,nNodesY,nNodesZ)
        i79  = jkl_to_global(j+4, k  , l+3,nNodesX,nNodesY,nNodesZ)
        i80  = jkl_to_global(j  , k+1, l+3,nNodesX,nNodesY,nNodesZ)
        i81  = jkl_to_global(j+1, k+1, l+3,nNodesX,nNodesY,nNodesZ)
        i82  = jkl_to_global(j+2, k+1, l+3,nNodesX,nNodesY,nNodesZ)
        i83  = jkl_to_global(j+3, k+1, l+3,nNodesX,nNodesY,nNodesZ)
        i84  = jkl_to_global(j+4, k+1, l+3,nNodesX,nNodesY,nNodesZ)
        i85  = jkl_to_global(j  , k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i86  = jkl_to_global(j+1, k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i87  = jkl_to_global(j+2, k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i88  = jkl_to_global(j+3, k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i89  = jkl_to_global(j+4, k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i90  = jkl_to_global(j  , k+3, l+3,nNodesX,nNodesY,nNodesZ)
        i91  = jkl_to_global(j+1, k+3, l+3,nNodesX,nNodesY,nNodesZ)
        i92  = jkl_to_global(j+2, k+3, l+3,nNodesX,nNodesY,nNodesZ)
        i93  = jkl_to_global(j+3, k+3, l+3,nNodesX,nNodesY,nNodesZ)
        i94  = jkl_to_global(j+4, k+3, l+3,nNodesX,nNodesY,nNodesZ)
        i95  = jkl_to_global(j  , k+4, l+3,nNodesX,nNodesY,nNodesZ)
        i96  = jkl_to_global(j+1, k+4, l+3,nNodesX,nNodesY,nNodesZ)
        i97  = jkl_to_global(j+2, k+4, l+3,nNodesX,nNodesY,nNodesZ)
        i98  = jkl_to_global(j+3, k+4, l+3,nNodesX,nNodesY,nNodesZ)
        i99  = jkl_to_global(j+4, k+4, l+3,nNodesX,nNodesY,nNodesZ)

        i100 = jkl_to_global(j  , k  , l+4,nNodesX,nNodesY,nNodesZ)
        i101 = jkl_to_global(j+1, k  , l+4,nNodesX,nNodesY,nNodesZ)
        i102 = jkl_to_global(j+2, k  , l+4,nNodesX,nNodesY,nNodesZ)
        i103 = jkl_to_global(j+3, k  , l+4,nNodesX,nNodesY,nNodesZ)
        i104 = jkl_to_global(j+4, k  , l+4,nNodesX,nNodesY,nNodesZ)
        i105 = jkl_to_global(j  , k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i106 = jkl_to_global(j+1, k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i107 = jkl_to_global(j+2, k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i108 = jkl_to_global(j+3, k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i109 = jkl_to_global(j+4, k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i110 = jkl_to_global(j  , k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i111 = jkl_to_global(j+1, k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i112 = jkl_to_global(j+2, k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i113 = jkl_to_global(j+3, k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i114 = jkl_to_global(j+4, k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i115 = jkl_to_global(j  , k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i116 = jkl_to_global(j+1, k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i117 = jkl_to_global(j+2, k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i118 = jkl_to_global(j+3, k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i119 = jkl_to_global(j+4, k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i120 = jkl_to_global(j  , k+4, l+4,nNodesX,nNodesY,nNodesZ)
        i121 = jkl_to_global(j+1, k+4, l+4,nNodesX,nNodesY,nNodesZ)
        i122 = jkl_to_global(j+2, k+4, l+4,nNodesX,nNodesY,nNodesZ)
        i123 = jkl_to_global(j+3, k+4, l+4,nNodesX,nNodesY,nNodesZ)
        i124 = jkl_to_global(j+4, k+4, l+4,nNodesX,nNodesY,nNodesZ)

        # first cube

        ee = ee + 1
        eConn[ee,:] = [i0, i50, i2, i10, i25, i26, i1, i30, i6, i5]

        ee = ee + 1
        eConn[ee,:] = [i52, i2, i50, i62, i27, i26, i51, i32, i56, i57]

        ee = ee + 1
        eConn[ee,:] = [i12, i62, i10, i2, i37, i36, i11, i32, i6, i7]

        ee = ee + 1
        eConn[ee,:] = [i60, i10, i62, i50, i35, i36, i61, i30, i56, i55]

        ee = ee + 1
        eConn[ee,:] = [i2, i50, i62, i10, i26, i56, i32, i30, i36, i6 ]

        # second cube

        ee = ee + 1
        eConn[ee,:] = [i52, i54, i2 , i62, i53, i28, i27, i58, i32, i57]

        ee = ee + 1
        eConn[ee,:] = [i4 , i2 , i54, i14, i3 , i28, i29, i8 , i34, i9]

        ee = ee + 1
        eConn[ee,:] = [i12, i14, i62, i2 , i13, i38, i37, i8 , i32, i7]

        ee = ee + 1
        eConn[ee,:] = [i64, i62, i14, i54, i63, i38, i39, i58, i34, i59]

        ee = ee + 1
        eConn[ee,:] = [i2 , i62, i54, i14, i32, i58, i28, i38, i34, i8]

        # third cube
      
        ee = ee + 1
        eConn[ee,:] = [i52 , i50 , i102, i62 , i51 , i76 , i77 , i56 , i82 , i57]

        ee = ee + 1
        eConn[ee,:] = [i100, i102, i50 , i110, i101, i76 , i75 , i106, i80 , i105]

        ee = ee + 1
        eConn[ee,:] = [i60 , i62 , i110, i50 , i61 , i86 , i85 , i56 , i80 , i55]

        ee = ee + 1
        eConn[ee,:] = [i112, i110, i62 , i102, i111, i86 , i87 , i106, i82 , i107]

        ee = ee + 1
        eConn[ee,:] = [i110, i102, i50 , i62 , i106, i76 , i80 , i82 , i56 , i86]

        # fourth cube
      
        ee = ee + 1
        eConn[ee,:] = [i52 , i102, i54 , i62 , i77 , i78 , i53 , i82 , i58 , i57]

        ee = ee + 1
        eConn[ee,:] = [i104, i54 , i102, i114, i79 , i78 , i103, i84 , i108, i109]

        ee = ee + 1
        eConn[ee,:] = [i64 , i114, i62 , i54 , i89 , i88 , i63 , i84 , i58 , i59]

        ee = ee + 1
        eConn[ee,:] = [i112, i62 , i114, i102, i87 , i88 , i113, i82 , i108, i107]

        ee = ee + 1
        eConn[ee,:] = [i54 , i102, i114, i62 , i78 , i108, i84 , i82 , i88 , i58]

        ee = ee + 1
        eConn[ee,:] = [i12 , i10 , i62 , i22 , i11 , i36 , i37 , i16 , i42 , i17]

        ee = ee + 1
        eConn[ee,:] = [i60 , i62 , i10 , i70 , i61 , i36 , i35 , i66 , i40 , i65 ]

        ee = ee + 1
        eConn[ee,:] = [i20 , i22 , i70 , i10 , i21 , i46 , i45 , i16 , i40 , i15]

        ee = ee + 1
        eConn[ee,:] = [i72 , i70 , i22 , i62 , i71 , i46 , i47 , i66 , i42 , i67]

        ee = ee + 1
        eConn[ee,:] = [i10 , i62 , i22 , i70 , i36 , i42 , i16 , i66 , i46 , i40]

        ee = ee + 1
        eConn[ee,:] = [i12 , i62 , i14 , i22 , i37 , i38 , i13 , i42 , i18 , i17]

        ee = ee + 1
        eConn[ee,:] = [i64 , i14 , i62 , i74 , i39 , i38 , i63 , i44 , i68 , i69 ]

        ee = ee + 1
        eConn[ee,:] = [i72 , i22 , i74 , i62 , i47 , i48 , i73 , i42 , i68 , i67]

        ee = ee + 1
        eConn[ee,:] = [i24 , i74 , i22 , i14 , i49 , i48 , i23 , i44 , i18 , i19]

        ee = ee + 1
        eConn[ee,:] = [i14 , i62 , i74 , i22 , i38 , i68 , i44 , i42 , i48 , i18]

        ee = ee + 1
        eConn[ee,:] = [i60 , i110, i62 , i70 , i85 , i86 , i61 , i90 , i66 , i65]

        ee = ee + 1
        eConn[ee,:] = [i112, i62 , i110, i122, i87 , i86 , i111, i92 , i116, i117]

        ee = ee + 1
        eConn[ee,:] = [i72 , i122, i70 , i62 , i97 , i96 , i71 , i92 , i66 , i67 ]

        ee = ee + 1
        eConn[ee,:] = [i120, i70 , i122, i110, i95 , i96 , i121, i90 , i116, i115]

        ee = ee + 1
        eConn[ee,:] = [i62 , i110, i122, i70 , i86 , i116, i92 , i90 , i96 , i66]

        ee = ee + 1
        eConn[ee,:] = [i64 , i62 , i114, i74 , i63 , i88 , i89 , i68 , i94 , i69 ]

        ee = ee + 1
        eConn[ee,:] = [i112, i114, i62 , i122, i113, i88 , i87 , i118, i92 , i117]

        ee = ee + 1
        eConn[ee,:] = [i72 , i74 , i122, i62 , i73 , i98 , i97 , i68 , i92 , i67]

        ee = ee + 1
        eConn[ee,:] = [i124, i122, i74 , i114, i123, i98 , i99 , i118, i94 , i119]

        ee = ee + 1
        eConn[ee,:] = [i62 , i114, i74 , i122, i88 , i94 , i68 , i118, i98 , i92]

      end # j-loop
      end # k-loop
      end # l-loop
  
    else    #  mesh several tetrahedra into corners
      
      # preallocate storage
      nElems = (nNodesX-3)*(nNodesY-3)*(nNodesZ-3)/64 * 40
      eConn = zeros(nElems,10)
      
      ee = 0
      for l=0:4:nNodesZ-3
      for k=0:4:nNodesY-3
      for j=0:4:nNodesX-3

        i1   = jkl_to_global(j  , k  , l  ,nNodesX,nNodesY,nNodesZ)
        i2   = jkl_to_global(j+1, k  , l  ,nNodesX,nNodesY,nNodesZ)
        i3   = jkl_to_global(j+2, k  , l  ,nNodesX,nNodesY,nNodesZ)
        i4   = jkl_to_global(j+3, k  , l  ,nNodesX,nNodesY,nNodesZ)
        i5   = jkl_to_global(j+4, k  , l  ,nNodesX,nNodesY,nNodesZ)
        i6   = jkl_to_global(j  , k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i7   = jkl_to_global(j+1, k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i8   = jkl_to_global(j+2, k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i9   = jkl_to_global(j+3, k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i10  = jkl_to_global(j+4, k+1, l  ,nNodesX,nNodesY,nNodesZ)
        i11  = jkl_to_global(j  , k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i12  = jkl_to_global(j+1, k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i13  = jkl_to_global(j+2, k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i14  = jkl_to_global(j+3, k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i15  = jkl_to_global(j+4, k+2, l  ,nNodesX,nNodesY,nNodesZ)
        i16  = jkl_to_global(j  , k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i17  = jkl_to_global(j+1, k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i18  = jkl_to_global(j+2, k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i19  = jkl_to_global(j+3, k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i20  = jkl_to_global(j+4, k+3, l  ,nNodesX,nNodesY,nNodesZ)
        i21  = jkl_to_global(j  , k+4, l  ,nNodesX,nNodesY,nNodesZ)
        i22  = jkl_to_global(j+1, k+4, l  ,nNodesX,nNodesY,nNodesZ)
        i23  = jkl_to_global(j+2, k+4, l  ,nNodesX,nNodesY,nNodesZ)
        i24  = jkl_to_global(j+3, k+4, l  ,nNodesX,nNodesY,nNodesZ)
        i25  = jkl_to_global(j+4, k+4, l  ,nNodesX,nNodesY,nNodesZ)

        i26  = jkl_to_global(j  , k  , l+1,nNodesX,nNodesY,nNodesZ)
        i27  = jkl_to_global(j+1, k  , l+1,nNodesX,nNodesY,nNodesZ)
        i28  = jkl_to_global(j+2, k  , l+1,nNodesX,nNodesY,nNodesZ)
        i29  = jkl_to_global(j+3, k  , l+1,nNodesX,nNodesY,nNodesZ)
        i30  = jkl_to_global(j+4, k  , l+1,nNodesX,nNodesY,nNodesZ)
        i31  = jkl_to_global(j  , k+1, l+1,nNodesX,nNodesY,nNodesZ)
#       i31.5= jkl_to_global(j+1, k+1, l+1,nNodesX,nNodesY,nNodesZ)
        i32  = jkl_to_global(j+2, k+1, l+1,nNodesX,nNodesY,nNodesZ)
#       i32.5= jkl_to_global(j+3, k+1, l+1,nNodesX,nNodesY,nNodesZ)
        i33  = jkl_to_global(j+4, k+1, l+1,nNodesX,nNodesY,nNodesZ)
        i34  = jkl_to_global(j  , k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i35  = jkl_to_global(j+1, k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i36  = jkl_to_global(j+2, k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i37  = jkl_to_global(j+3, k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i38  = jkl_to_global(j+4, k+2, l+1,nNodesX,nNodesY,nNodesZ)
        i39  = jkl_to_global(j  , k+3, l+1,nNodesX,nNodesY,nNodesZ)
#       i39.5= jkl_to_global(j+1, k+3, l+1,nNodesX,nNodesY,nNodesZ)
        i40  = jkl_to_global(j+2, k+3, l+1,nNodesX,nNodesY,nNodesZ)
#       i40.5= jkl_to_global(j+3, k+3, l+1,nNodesX,nNodesY,nNodesZ)
        i41  = jkl_to_global(j+4, k+3, l+1,nNodesX,nNodesY,nNodesZ)
        i42  = jkl_to_global(j  , k+4, l+1,nNodesX,nNodesY,nNodesZ)
        i43  = jkl_to_global(j+1, k+4, l+1,nNodesX,nNodesY,nNodesZ)
        i44  = jkl_to_global(j+2, k+4, l+1,nNodesX,nNodesY,nNodesZ)
        i45  = jkl_to_global(j+3, k+4, l+1,nNodesX,nNodesY,nNodesZ)
        i46  = jkl_to_global(j+4, k+4, l+1,nNodesX,nNodesY,nNodesZ)

        i47  = jkl_to_global(j  , k  , l+2,nNodesX,nNodesY,nNodesZ)
        i48  = jkl_to_global(j+1, k  , l+2,nNodesX,nNodesY,nNodesZ)
        i49  = jkl_to_global(j+2, k  , l+2,nNodesX,nNodesY,nNodesZ)
        i50  = jkl_to_global(j+3, k  , l+2,nNodesX,nNodesY,nNodesZ)
        i51  = jkl_to_global(j+4, k  , l+2,nNodesX,nNodesY,nNodesZ)
        i52  = jkl_to_global(j  , k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i53  = jkl_to_global(j+1, k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i54  = jkl_to_global(j+2, k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i55  = jkl_to_global(j+3, k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i56  = jkl_to_global(j+4, k+1, l+2,nNodesX,nNodesY,nNodesZ)
        i57  = jkl_to_global(j  , k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i58  = jkl_to_global(j+1, k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i59  = jkl_to_global(j+2, k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i60  = jkl_to_global(j+3, k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i61  = jkl_to_global(j+4, k+2, l+2,nNodesX,nNodesY,nNodesZ)
        i62  = jkl_to_global(j  , k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i63  = jkl_to_global(j+1, k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i64  = jkl_to_global(j+2, k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i65  = jkl_to_global(j+3, k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i66  = jkl_to_global(j+4, k+3, l+2,nNodesX,nNodesY,nNodesZ)
        i67  = jkl_to_global(j  , k+4, l+2,nNodesX,nNodesY,nNodesZ)
        i68  = jkl_to_global(j+1, k+4, l+2,nNodesX,nNodesY,nNodesZ)
        i69  = jkl_to_global(j+2, k+4, l+2,nNodesX,nNodesY,nNodesZ)
        i70  = jkl_to_global(j+3, k+4, l+2,nNodesX,nNodesY,nNodesZ)
        i71  = jkl_to_global(j+4, k+4, l+2,nNodesX,nNodesY,nNodesZ)

        i72  = jkl_to_global(j  , k  , l+3,nNodesX,nNodesY,nNodesZ)
        i73  = jkl_to_global(j+1, k  , l+3,nNodesX,nNodesY,nNodesZ)
        i74  = jkl_to_global(j+2, k  , l+3,nNodesX,nNodesY,nNodesZ)
        i75  = jkl_to_global(j+3, k  , l+3,nNodesX,nNodesY,nNodesZ)
        i76  = jkl_to_global(j+4, k  , l+3,nNodesX,nNodesY,nNodesZ)
        i77  = jkl_to_global(j  , k+1, l+3,nNodesX,nNodesY,nNodesZ)
#       i77.5= jkl_to_global(j+1, k+1, l+3,nNodesX,nNodesY,nNodesZ)
        i78  = jkl_to_global(j+2, k+1, l+3,nNodesX,nNodesY,nNodesZ)
#       i78.5= jkl_to_global(j+3, k+1, l+3,nNodesX,nNodesY,nNodesZ)
        i79  = jkl_to_global(j+4, k+1, l+3,nNodesX,nNodesY,nNodesZ)
        i80  = jkl_to_global(j  , k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i81  = jkl_to_global(j+1, k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i82  = jkl_to_global(j+2, k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i83  = jkl_to_global(j+3, k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i84  = jkl_to_global(j+4, k+2, l+3,nNodesX,nNodesY,nNodesZ)
        i85  = jkl_to_global(j  , k+3, l+3,nNodesX,nNodesY,nNodesZ)
#       i85.5= jkl_to_global(j+1, k+3, l+3,nNodesX,nNodesY,nNodesZ)
        i86  = jkl_to_global(j+2, k+3, l+3,nNodesX,nNodesY,nNodesZ)
#       i86.5= jkl_to_global(j+3, k+3, l+3,nNodesX,nNodesY,nNodesZ)
        i87  = jkl_to_global(j+4, k+3, l+3,nNodesX,nNodesY,nNodesZ)
        i88  = jkl_to_global(j  , k+4, l+3,nNodesX,nNodesY,nNodesZ)
        i89  = jkl_to_global(j+1, k+4, l+3,nNodesX,nNodesY,nNodesZ)
        i90  = jkl_to_global(j+2, k+4, l+3,nNodesX,nNodesY,nNodesZ)
        i91  = jkl_to_global(j+3, k+4, l+3,nNodesX,nNodesY,nNodesZ)
        i92  = jkl_to_global(j+4, k+4, l+3,nNodesX,nNodesY,nNodesZ)

        i93  = jkl_to_global(j  , k  , l+4,nNodesX,nNodesY,nNodesZ)
        i94  = jkl_to_global(j+1, k  , l+4,nNodesX,nNodesY,nNodesZ)
        i95  = jkl_to_global(j+2, k  , l+4,nNodesX,nNodesY,nNodesZ)
        i96  = jkl_to_global(j+3, k  , l+4,nNodesX,nNodesY,nNodesZ)
        i97  = jkl_to_global(j+4, k  , l+4,nNodesX,nNodesY,nNodesZ)
        i98  = jkl_to_global(j  , k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i99  = jkl_to_global(j+1, k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i100 = jkl_to_global(j+2, k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i101 = jkl_to_global(j+3, k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i102 = jkl_to_global(j+4, k+1, l+4,nNodesX,nNodesY,nNodesZ)
        i103 = jkl_to_global(j  , k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i104 = jkl_to_global(j+1, k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i105 = jkl_to_global(j+2, k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i106 = jkl_to_global(j+3, k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i107 = jkl_to_global(j+4, k+2, l+4,nNodesX,nNodesY,nNodesZ)
        i108 = jkl_to_global(j  , k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i109 = jkl_to_global(j+1, k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i110 = jkl_to_global(j+2, k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i111 = jkl_to_global(j+3, k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i112 = jkl_to_global(j+4, k+3, l+4,nNodesX,nNodesY,nNodesZ)
        i113 = jkl_to_global(j  , k+4, l+4,nNodesX,nNodesY,nNodesZ)
        i114 = jkl_to_global(j+1, k+4, l+4,nNodesX,nNodesY,nNodesZ)
        i115 = jkl_to_global(j+2, k+4, l+4,nNodesX,nNodesY,nNodesZ)
        i116 = jkl_to_global(j+3, k+4, l+4,nNodesX,nNodesY,nNodesZ)
        i117 = jkl_to_global(j+4, k+4, l+4,nNodesX,nNodesY,nNodesZ)

        # first cube
        ee = ee + 1
        eConn[ee,:] = [i47 , i49 , i1  , i57 , i48 , i27 , i26 , i53 , i31 , i52 ]

        ee = ee + 1
        eConn[ee,:] = [i3  , i1  , i49 , i13 , i2  , i27 , i28 , i7  , i32 , i8  ]

        ee = ee + 1
        eConn[ee,:] = [i11 , i13 , i57 , i1  , i12 , i35 , i34 , i7  , i31 , i6  ]

        ee = ee + 1
        eConn[ee,:] = [i59 , i57 , i13 , i49 , i58 , i35 , i36 , i53 , i32 , i54 ]

        ee = ee + 1
        eConn[ee,:] = [i1  , i57 , i49 , i13 , i31 , i53 , i27 , i35 , i32 , i7  ]

        # second cube
        ee = ee + 1
        eConn[ee,:] = [i3  , i49 , i5  , i13 , i28 , i29 , i4  , i32 , i9  , i8  ]

        ee = ee + 1
        eConn[ee,:] = [i51 , i5  , i49 , i61 , i30 , i29 , i50 , i33 , i55 , i56 ]

        ee = ee + 1
        eConn[ee,:] = [i15 , i61 , i13 , i5  , i38 , i37 , i14 , i33 , i9  , i10  ]

        ee = ee + 1
        eConn[ee,:] = [i59 , i13 , i61 , i49 , i36 , i37 , i60 , i32 , i55 , i54 ]

        ee = ee + 1
        eConn[ee,:] = [i13 , i49 , i5  , i61 , i32 , i29 , i9  , i55 , i33 , i37 ]

        # third cube
        ee = ee + 1
        eConn[ee,:] = [i11 , i57 , i13 , i21 , i34 , i35 , i12 , i39 , i17 , i16 ]

        ee = ee + 1
        eConn[ee,:] = [i59 , i13 , i57 , i69 , i36 , i35 , i58 , i40 , i63 , i64 ]

        ee = ee + 1
        eConn[ee,:] = [i23 , i69 , i21 , i13 , i44 , i43 , i22 , i40 , i17 , i18 ]

        ee = ee + 1
        eConn[ee,:] = [i67 , i21 , i69 , i57 , i42 , i43 , i68 , i39 , i63 , i62 ]

        ee = ee + 1
        eConn[ee,:] = [i13 , i69 , i21 , i57 , i40 , i43 , i17 , i63 , i39 , i35 ]

        # fourth cube
        ee = ee + 1
        eConn[ee,:] = [i15 , i13 , i61 , i25 , i14 , i37 , i38 , i19 , i41 , i20 ]

        ee = ee + 1
        eConn[ee,:] = [i59 , i61 , i13 , i69 , i60 , i37 , i36 , i65 , i40 , i64 ]

        ee = ee + 1
        eConn[ee,:] = [i23 , i25 , i69 , i13 , i24 , i45 , i44 , i19 , i40 , i18 ]

        ee = ee + 1
        eConn[ee,:] = [i71 , i69 , i25 , i61 , i70 , i45 , i46 , i65 , i41 , i66 ]

        ee = ee + 1
        eConn[ee,:] = [i13 , i69 , i61 , i25 , i40 , i65 , i37 , i45 , i41 , i19 ]

        # fifth cube
        ee = ee + 1
        eConn[ee,:] = [i47 , i93 , i49 , i57 , i72 , i73 , i48 , i77 , i53 , i52 ]

        ee = ee + 1
        eConn[ee,:] = [i95 , i49 , i93 , i105, i74 , i73 , i94 , i78 , i99 , i100]

        ee = ee + 1
        eConn[ee,:] = [i59 , i105 , i57 , i49 , i82 , i81 , i58 , i78 , i53 , i54 ]

        ee = ee + 1
        eConn[ee,:] = [i103, i57 , i105, i93 , i80 , i81 , i104, i77 , i99 , i98 ]

        ee = ee + 1
        eConn[ee,:] = [i49 , i57 , i93 , i105, i53 , i77 , i73 , i81 , i99 , i78 ]

        # sixth cube
        ee = ee + 1
        eConn[ee,:] = [i51 , i49 , i97 , i61 , i50 , i75 , i76 , i55 , i79 , i56 ]

        ee = ee + 1
        eConn[ee,:] = [i95 , i97 , i49 , i105, i96 , i75 , i74 , i101, i78 , i100]

        ee = ee + 1
        eConn[ee,:] = [i59 , i61 , i105, i49 , i60 , i83 , i82 , i55 , i78 , i54 ]

        ee = ee + 1
        eConn[ee,:] = [i107, i105, i61 , i97 , i106, i83 , i84 , i101, i79 , i102]

        ee = ee + 1
        eConn[ee,:] = [i49 , i61 , i105, i97 , i55 , i83 , i78 , i79 , i101, i75 ]

        # seventh cube
        ee = ee + 1
        eConn[ee,:] = [i59 , i57 , i105, i69 , i58 , i81 , i82 , i63 , i86 , i64 ]

        ee = ee + 1
        eConn[ee,:] = [i103, i105, i57 , i113, i104, i81 , i80 , i109, i85 , i108]

        ee = ee + 1
        eConn[ee,:] = [i67 , i69 , i113, i57 , i68 , i89 , i88 , i63 , i85 , i62 ]

        ee = ee + 1
        eConn[ee,:] = [i115, i113, i69 , i105, i114, i89 , i90 , i109, i86 , i110]
  
        ee = ee + 1
        eConn[ee,:] = [i57 , i113, i105, i69 , i85 , i109, i81 , i89 , i86 , i63 ]

        # eighth cube
        ee = ee + 1
        eConn[ee,:] = [i59 , i105, i61 , i69 , i82 , i83 , i60 , i86 , i65 , i64 ]

        ee = ee + 1
        eConn[ee,:] = [i107, i61 , i105, i117, i84 , i83 , i106, i87 , i111, i112]

        ee = ee + 1
        eConn[ee,:] = [i71 , i117, i69 , i61 , i92 , i91 , i70 , i87 , i65 , i66 ]

        ee = ee + 1
        eConn[ee,:] = [i115, i69 , i117, i105, i90 , i91 , i116, i86 , i111, i110]

        ee = ee + 1
        eConn[ee,:] = [i61 , i105, i117, i69 , i83 , i111, i87 , i86 , i91 , i65 ]

      end # j-loop
      end # k-loop
      end # l-loop
    end
  
  else
    error("TetMesh_cubeMesh: Orders other than 1 or 2 are not implemented \n")
  end

  return x, eConn
end


function jkl_to_global(j,k,l,nNodesX,nNodesY,~)

  # use logic from c-code (j,k,l indices start with 0)

  l_odd = mod(l,2)

  be = (l+l_odd)/2    # number of preceeding even numbered l planes
  bo = (l-l_odd)/2    # "                    odd  "               "

  basel = be*( nNodesX*nNodesY ) + bo*( nNodesX*(nNodesY+1)/2 + (nNodesX+1)/2*(nNodesY-1)/2 )

  k_odd = mod(k,2)

  if l_odd
    base = basel + (k+k_odd)/2 *nNodesX + (k-k_odd)/2 *(nNodesX+1)/2
  else
    base = basel + k*nNodesX
  end

  if k_odd && l_odd     # skip center of cubes
    i = j/2 + base
  else
    i = j + base
  end

  i = i+1  # 1-based indexing

  return i
end
