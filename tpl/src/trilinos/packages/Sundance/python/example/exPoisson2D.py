#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance


import math
from PySundance import *

class LeftPointPredicate :
  def evalOp(self, x, y) :
      return (math.fabs(x) < 1.0e-10);

class RightPointPredicate :
  def evalOp(self, x, y) :
      return (math.fabs(x-1.0) < 1.0e-10);

class BottomPointPredicate :
  def evalOp(self, x, y) :
      return (math.fabs(y) < 1.0e-10);

class TopPointPredicate :
  def evalOp(self, x, y) :
      return (math.fabs(y-2.0) < 1.0e-10);

def main():

  vecType = EpetraVectorType()
  npx = 1
  npy = getNProc()
  ny = 32
  nx = 32
  mesher  = PartitionedRectangleMesher(0.0, 1.0, nx, npx,
                                       0.0, 2.0, ny/npy, npy);
  mesh = mesher.getMesh();

  check = mesh.checkConsistency('meshCheck')


  if check==0 :
    print 'INCONSISTENT MESH'
  basis = Lagrange(2)

  u = UnknownFunction(basis, "u");
  v = TestFunction(basis, "v");
  x = CoordExpr(0);
  y = CoordExpr(1);
  dx = Derivative(0);
  dy = Derivative(1);
  grad = List(dx, dy);
  
  quad2 = GaussianQuadrature(2)
  quad4 = GaussianQuadrature(4)
  
  interior = MaximalCellFilter()
  edges = DimensionalCellFilter(1)
  left = edges.subset(LeftPointPredicate())
  right = edges.subset(RightPointPredicate())
  top = edges.subset(TopPointPredicate())
  bottom = edges.subset(BottomPointPredicate())

  one = 1.0
  eqn = Integral(interior, (grad*v)*(grad*u) + one*v, quad2)\
        + Integral(top, -v/3.0, quad2)\
        + Integral(right, -v*(1.5 + (1.0/3.0)*y - u), quad4)
  
  bc = EssentialBC(bottom, v*(u-0.5*x*x), quad4)
  
  prob = LinearProblem(mesh, eqn, bc, v, u, vecType)

  solver = readSolver(searchForFile("aztec-ml.xml"));

  soln = prob.solve(solver)


  exactSoln = 0.5*x*x + (1.0/3.0)*y;

  writer = VTKWriter("Poisson2D");
  writer.addMesh(mesh)
  writer.addField("u0", soln)
  writer.write()
  

  diff = (soln - exactSoln)**2.0
  diffDeriv = (dx*(soln - exactSoln))**2.0
  diffFlux = (dy*(soln - exactSoln))**2.0
  
  error = math.sqrt(diff.integral(interior, mesh, quad4))
  derivError = math.sqrt(diffDeriv.integral(interior, mesh, quad4))
  fluxError = math.sqrt(diffFlux.integral(top, mesh, quad4))
  print "error = " , error
  print "deriv error = " , derivError
  print "flux error = " , fluxError

  error = max(error, derivError)
  error = max(error, fluxError)

  print "max err = ", error

  tol = 1.0e-8
  passFailTest(error, tol)



  
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
