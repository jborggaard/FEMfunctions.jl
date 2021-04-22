# FEMfunctions

[![Build Status](https://travis-ci.com/jborggaard/FEMfunctions.jl.svg?branch=master)](https://travis-ci.com/jborggaard/FEMfunctions.jl)
[![Coverage](https://codecov.io/gh/jborggaard/FEMfunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jborggaard/FEMfunctions.jl)


This collection of *julia* functions are used for _oned_, _twod_, or _threed_ 
PDE solvers based on the finite element method (FEM).  A list of available 
functions with brief descriptions, followed by an explanation of calling  arguments for each functions follows below (grouped by dimension).

To install:
```julia
  julia> using Pkg
  julia> Pkg.add(url="https://github.com/jborggaard/FEMfunctions.jl")
````

Copyright 2013 Jeff Borggaard, 
Department of Mathematics,
Interdisciplinary Center for Applied Mathematics,
Virginia Tech

This file is part of the FEMfunctions package.

This library of functions is free software: you can redistribute it and/or 
modify it under the terms of the MIT License.

FEMfunctions is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.

---
## oned functions ##


```julia
  A = onedBilinear(kernel,phi,test,quadratureWeights);
````
- `onedBilinear`: Integrates bilinear form over an element
  by computing \int{ kernel . phi . test }
  
   Inputs | Description
   ------ | -----------
   kernel | a kernel function evaluated at quadrature points
   phi    | matrix of basis functions evaluated at quadrature points
   test   | (test) basis functions evaluated at quadrature points
   quadratureWeights | ... scaled by Jacobian evaluated at quadrature points
   
   Outputs | Description
   ------- | -----------
   A       | \int{ kernel . phi . test } (nTest x nPhi matrix)
          
  
- `onedLinForm`: Integrates linear forms: functionals over an element

  ```julia
    f = onedLinForm(kernel,test,quadratureWeights)
  ````
  routine to compute \int{ kernel*test }

- `onedQuadrature`: Provides a few Gaussian quadrature rules

  ```julia
    r,w = onedQuadrature(rule)
  ````
  calculate quadrature integration points on (-1,1)
  
- `onedGalerkinProjection`: Computes Galerkin projection of a given function onto a finite element mesh (useful for setting initial conditions)

  ```julia
    p = onedGalerkinProjection(x,eConn,function)
  ````
  routine to perform Galerkin projection of the incoming function.
  
- `onedMesh`: Provides a one-dimensional mesh given a mesh density

  ```julia
  x,eConn,indexU,indexD = oned_mesh(xb,eConnb,rho)
  ````
   Inputs | Description
   ------ | -----------
   xb     | nodes for a background mesh
   eConnb | element connectivity for background mesh
   rho    | a vector of element densities per element
   
   This is frequently used to generate a uniform mesh with equispaced points, e.g.
   ```julia
      [x,eConn,indexU,indexD] = onedMesh([0.0 1.0],[1 2],[20])
   ````

- `onedProjectDerivative`: Projects the derivative of a finite element solution onto a continuous finite element space

  ```julia
    d, elementError = onedProjectDerivative(x,eConn,u)
  ````

- `onedShape`: Evaluates finite element shape functions in an element at the quadrature points

  ```julia
    [xg,wg,phi,p_x,p_xx] = onedShape(xLocal,r,w)
  ````
   Inputs | Description
   ------ | -----------
   xLocal | nodes for one element
   r,w    | output from onedQuadrature
   
   Outputs | Description
   ------- | -----------
   xg,wg   | points and weights for current element
   phi     | values of finite element shape functions at xg
   p_x     | derivative of finite element shape function at xg
   p_xx    | second derivative (if x has 3 or more nodes)
   

- `onedShapeHermite`: Same as `onedShape`, except computing Hermite finite element shape functions within one element at quadrature points

  ```julia
    [xg,wg,phi0,phi1,p0_x,p1_x,p0_xx,p1_xx] = onedShapeHermite(xLocal,r,w)
  ````

- `onedShapeIso`: Same as `onedShape`, but allows for isoparametric
elements (interior nodes are not uniformly spaced)

  ```julia
    [xg,wg,phi,p_x,p_xx] = oned_shapeiso(xLocal,r,w)
  ````

---
## twod functions ##

These are functions for P# (triangular) or Q# (rectangular) elements.

- `twodBilinear`: Integrates bilinear form over an element

  ```julia
    A = twodBilinear( kernel, phi, test, wg )
  ````

- `twodLinForm`: Integrates functionals over an element

  ```julia
    F = twodLinForm(kernel, test, wg)
  ````

- `twodQuadratureQ`: Provides a few quadrature rules for Q# elements

  ```julia
    [r,s,w] = twodQuadratureQ(rule)
  ````

### Specific functions for triangular meshes P# ###
- `twodGetFaces`: Uses element connectivity to find boundary faces

  ```julia
    elements,faces = twodGetFaces( eConn, index )
  ````
   Inputs | Description
   ------ | -----------
   eConn  | element connectivity
   index  | searches over a subset of elements (optional)

- `twodMesh`: Builds a triangular mesh on a rectangular domain

  ```julia
    x,eConn,indexBoundary = twodMesh(xl, xr, yl, yr, 
                                     elementType,       
                                     nNodesx, nNodesy)
  ````

- `twodOrientate`: Reorders node numbers to 

  ```julia
    t = twodOrientate(p, t)
  ````

- `twodQuadrature`: Provides a few quadrature rules for P# elements

  ```julia
    r,s,w = twodQuadrature(rule)
  ````

- `twodShape`: Evaluates finite element shape functions in an element at the quadrature points

  ```julia
    xg, wg, phi, p_x, p_y = twodShape(xLocal, r, s, w)
  ````

- `twodShapeIso`: Same as twoShape except allows isoparametric elements

  ```julia
    xg, wg, phi, p_x, p_y = twodShapeIso(xLocal, r, s, w)
  ````

  --
- `TriMesh_ElementAdjacency`: Returns the adjacency graph of a mesh

  ```julia
    eAdjacency = TriMesh_ElementAdjacency(eConn)   
  ````

- `TriMesh_Interpolate`:

  ```julia
    f,elementList = TriMesh_Interpolate( x, eConn, neighbors, f_node, xy_points)
  ````

- `TriMesh_ProjectDerivatives`: use patches to project derivatives of finite element functions onto a continuous finite element function

  ```julia
    d1_p, d2_p, e_error, node = twod_ProjectDerivatives(x, eConn, u)
  ````

- `TriMesh_PromoteL2Q`: promotes a linear finite element representation to a quadratic one

  ```julia
    uQuadratic = TriMesh_PromoteL2Q(uLinear, eConn)
  ````
  
- `TriMesh_RestrictQ2L`: demotes a quadratic finite element representation to a linear one

  ```julia
    x, eConn, uLinear = TriMesh_RestrictQ2L( x, eConn, uQuadratic )
  ````

- `TriMesh_Search`: searches for the element containing a specified point

  ```julia
    x, iso2, iso3 = TriMesh_Search(x, eConn, eAdjacency, point)
  ````
  
### Specific functions for rectangular meshes Q# ###
- `twodGridQ`: Creates a uniform grid of nodes marking the boundary nodes

  ```julia
	xy,boundaryNodes = twodGridQ(jmax, kmax, jdim, kdim)
  ````

- `twodMeshQ`: Creates a uniform 2D quadrateral mesh 

  ```julia
    x, eConn, bNodes = twodMeshQ( order, nx, ny, domain )
  ````

- `twodShapeQ`:

  ```julia
    xg, wg, phi, p_x, p_y = twodShapeQ(x, r, s, w)
  ````

### Functions to save FEM solutions for visualization ###
- `saveFEMasGmsh`:

  ```julia
    saveFEMasGmsh(filename, x, eConn, u, title)
  ````

- `saveFEMasVTK`:

  ```julia
    saveFEMasVTK(filename, 
                 x, eConn, 
                 scalarLabels, scalars, 
                 vectorLabels, vectors)
  ````

---
## threed functions ##

- `threedBilinear`: Integrates bilinear form over an element

  ```julia
    A = threedBilinear( kernel, phi, test, wg )
  ````

- `threedLinForm`: Integrates functionals over an element

  ```julia
    f = threedLinForm( kernel, test, wg )
  ````

- `threedQuadrature`: Provides a few quadrature rules

  ```julia
    r,s,t,w = threedQuadrature(rule)
  ````

- `threedMesh`: Generates a tetrahedral mesh for a cubic domain

  ```julia
    x, eConn = threed_mesh(order, xl,yl,zl, jmax,kmax,lmax)
  ````

- `threedOrientate`:

  ```julia
    eConnOrientated = threed_orientate( x,eConn )
  ````

- `threedProjectDerivatives`:

  ```julia
    d1_p, d2_p, d3_p, e_error, node] = threedProjectDerivatives(x,e_conn,u,nQuadrature)
  ````

- `threedRestrictQ2L`:

  ```julia
    x, eConn, uLinear = threedRestrictQ2L( x, eConn, uQuadratic )
  ````
  Some visualization tools require linear tetrahedral elements.

- `threedShape`:

  ```julia
    xg,wg,phi,p_x,p_y,p_z = threedShape(xLocal, r, s, t, w)
  ````

- `threedShapeQ`:

  ```julia
    xg,wg,phi,p_x,p_y,p_z = threedShapeQ(xLocal, r, s, t, w)
  ````

- `threedShapeIso`:

  ```julia
    xg,wg,phi,p_x,p_y,p_z = threedShapeIso(xLocal, r, s, t, w)
  ````



Author
--------
Jeff Borggaard, Interdisciplinary Center for Applied Mathematics, Virginia Tech
jborggaard@vt.edu

License
--------
These files are provided under the MIT License.

* Maintained on GitHub:  [jborggaard.github.com/fem_functions](jborggaard.github.com/FEMfunctions.jl)
