function onedElliptic(nElements=40,order=2)
#  Solves a linear elliptic PDE in 1D 

  f(x)=x           # right-hand side function definition.

  x,eConn,iU,iC = onedMesh( [0.0,1.0], order, nElements )
  nNodes = size(x,1)

  ide = [1:nNodes;]     # index into unknowns  (set b.c. later in 1D)

  r,w = onedQuadratureRule(5)
  o   = ones(size(r))

  nElDOF = order+1
  nEntries = nElements*nElDOF^2
  II  = zeros(Int64,nEntries) #Array{Int32,1}(undef,nEntries)     # zeros(Int32,nEntries)
  JJ  = zeros(Int64,nEntries) #Array{Int32,1}(undef,nEntries)     # zeros(Int32,nEntries)
  AA  = zeros(Float64,nEntries) #Array{Float64,1}(undef,nEntries)   # zeros(Float64,nEntries)
  b   = zeros(Float64,nNodes) #Array{Float64,1}(undef,nNodes)     # zeros(Float64,nNodes)

  index = 0
  for k=1:nElements
    nLocal = eConn[k,:]
    xLocal = view(x,nLocal)

    xg,wg,phi,p_x,p_xx = onedShape( xLocal, r, w )
    fg     = f(xg)             # forcing function evaluated at quadrature points

    ALoc   = onedBilinear(  o, p_x, p_x, wg )
    bLoc   = onedLinForm(  fg,      phi, wg )

    lDOF = ide[nLocal]
    for nt = 1:nElDOF
      for nu = 1:nElDOF
        index = index + 1
        II[index]   = lDOF[nu]
        JJ[index]   = lDOF[nt]
        AA[index]   = ALoc[nt,nu]
      end
      b[lDOF[nt]] = b[lDOF[nt]] + bLoc[nt]
    end
  end

  A = sparse(II,JJ,AA)

  u = A[iU,iU]\b[iU]

  #  Calculate the exact solution for unit test
  xi = x[iU];

  uExact = -(xi.^3-xi)/6

  return norm(u-uExact)
# using Gaston
# figure(1)
# Gaston.set_terminal("x11");
# Gaston.closeall();
# c = Gaston.CurveConf();
# c.color = "blue";
# c.legend = "FEM";
# Gaston.addcoords(xi,u,c);

# figure(2)
# Gaston.set_terminal("x11");
# c = Gaston.CurveConf();
# c.color = "red";
# c.legend = "Exact";
# c.marker = "fdmd";
# c.plotstyle = "linespoints";
# Gaston.addcoords(xi,-(xi.^3-xi)/6,c);
# Gaston.llplot();
  
# set_filename("test.png")
# printfigure("png")


#  using PyPlot
#  plot(xi,u,"r",label="FEM solution")
#  plot(xi,-(xi.^3-xi)/6,"go",label="Exact at nodes")
#  legend()
end
