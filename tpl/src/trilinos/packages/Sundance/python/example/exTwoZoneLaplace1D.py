#!/usr/bin/python -u

import setpath
import PySundance


import math
from PySundance import *
from aztecSolver import solverParams

class LeftPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt) < 1.0e-10);

class RightPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt-1.0) < 1.0e-10);

class LeftZonePredicate :
  def evalOp(self, pt) :
    rtn = (pt <= 0.5)
    return rtn;

class RightZonePredicate :
  def evalOp(self, pt) :
    rtn = (pt >= 0.5)
    return rtn;



def main():

  vecType = EpetraVectorType()
  nProc = getNProc()
  mesher  = PartitionedLineMesher(0.0, 1.0, 10*nProc);
  mesh = mesher.getMesh();
  basis = Lagrange(2)

  u = UnknownFunction(basis, "u");
  v = TestFunction(basis, "v");
  x = CoordExpr(0);
  dx = Derivative(0);

  quad = GaussianQuadrature(2)
  
  interior = MaximalCellFilter()
  points = DimensionalCellFilter(0)
  leftPt = points.subset(LeftPointPredicate())
  rightPt = points.subset(RightPointPredicate())
  leftZone = interior.subset(LeftZonePredicate())
  rightZone = interior.subset(RightZonePredicate())

  kappa1 = 1.0
  kappa2 = 0.5
  
  bc = EssentialBC(leftPt, v*u, quad) + EssentialBC(rightPt, v*(u-1.0), quad)
  eqn = Integral(leftZone, kappa1*(dx*v)*(dx*u), quad) \
        + Integral(rightZone, kappa2*(dx*v)*(dx*u), quad)
  
  prob = LinearProblem(mesh, eqn, bc, v, u, vecType)

  solver = readSolver(searchForFile("aztec.xml"));
  
  soln = prob.solve(solver)

  a1 = 2.0/(1 + kappa1/kappa2)
  a2 = kappa1*a1/kappa2

  exactLeft = a1*x
  exactRight = a2*(x-0.5) + a1/2.0
  diffLeft = 0.5*(soln-exactLeft)**2
  diffRight = 0.5*(soln-exactRight)**2

  error2 = diffLeft.integral(leftZone, mesh, quad) \
           + diffRight.integral(rightZone, mesh, quad)
  error = math.sqrt(error2)
  
  print "error = " , error


  tol = 1.0e-11
  passFailTest(error, tol)
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
