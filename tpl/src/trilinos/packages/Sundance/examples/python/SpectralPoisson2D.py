#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance

import math
from PySundance import *


# define solver parameters:
amesosSolverParams = ParameterList({"Linear Solver" : 
                                    {"Type" : "Amesos",
                                     "Kernel" : "Umfpack"
                                     }
                                    })

def main():
  pi = 4.0*math.atan(1.0)

  vecType = EpetraVectorType()
  npx = 1
  npy = getNProc()
  ny = 96
  nx = 96
  mesher  = PartitionedRectangleMesher(0.0, 1.0, nx, npx,
                                       0.0, 1.0, ny/npy, npy);
  mesh = mesher.getMesh();

  basis = Lagrange(1)
  spBasis = HermiteSpectralBasis(1, 2)

  u = UnknownFunction(basis, spBasis, "u");
  v = TestFunction(basis, spBasis, "v");
  x = CoordExpr(0);
  y = CoordExpr(1);
  dx = Derivative(0);
  dy = Derivative(1);
  grad = List(dx, dy);
  
  quad2 = GaussianQuadrature(2)
  quad4 = GaussianQuadrature(4)
  
  interior = MaximalCellFilter()
  bdry = BoundaryCellFilter()

  kappa = SpectralExpr(spBasis, List(Expr(1.0), 0.2*sin(pi*x)*sin(pi*y),
                                     0.1*sin(2.0*pi*x)*sin(2.0*pi*y)))
  
  eqn = Integral(interior, (grad*v)*((grad*u)*kappa) + v, quad4)

  print 'eqn = ', eqn
  
  bc = EssentialBC(bdry, v*u, quad4)
  
  prob = LinearProblem(mesh, eqn, bc, v, u, vecType)

  solver = buildSolver(amesosSolverParams)

  print 'solving...'

  soln = prob.solve(solver)

  writer = VTKWriter("Poisson2D");
  writer.addMesh(mesh)
  for i in range(spBasis.nterms()):
      writer.addField("u%d" % i, soln[i])
  writer.write()
  

  
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
