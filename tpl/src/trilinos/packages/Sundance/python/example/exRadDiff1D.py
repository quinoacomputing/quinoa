#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance


import math
from PySundance import *

class LeftPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt) < 1.0e-10);

class RightPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt-1.0) < 1.0e-10);

from noxSolver import solverParams
from noxSolver import solverDict

def main():

  vecType = EpetraVectorType()
  nx = 400
  mesher  = PartitionedLineMesher(0.0, 1.0, nx*getNProc());
  mesh = mesher.getMesh();
  basis = Lagrange(1)

  u = UnknownFunction(basis, "u");
  v = TestFunction(basis, "v");
  x = CoordExpr(0);
  dx = Derivative(0);
  
  quad = GaussianQuadrature(8)

  discSpace = DiscreteSpace(mesh, basis, vecType)
  projector = L2Projector(discSpace, 1.0 + x)
  u0 = projector.project()
  
  interior = MaximalCellFilter()
  points = DimensionalCellFilter(0)
  leftPt = points.subset(LeftPointPredicate())
  rightPt = points.subset(RightPointPredicate())
  
  bc = EssentialBC(leftPt, v*(u-1.0), quad) \
       + EssentialBC(rightPt, v*(u-2.0), quad)
  eqn = Integral(interior, u*u*u*(dx*v)*(dx*u), quad)

  prob = NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType)
  solver = NOXSolver(solverParams)

  prob.solve(solver)

  exactSoln = pow(15.0*x + 1.0, 0.25)

  diff = (u0 - exactSoln)**2

#  print 'soln=\n'
#  print u0.getVector()

  error = math.sqrt(diff.integral(interior, mesh, quad))
  print "error = " , error

  tol = 1.0e-5
  passFailTest(error, tol)
  
  
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
