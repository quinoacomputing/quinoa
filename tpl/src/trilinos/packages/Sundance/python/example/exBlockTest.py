#!/usr/bin/env python

import setpath
import PySundance


import math
from PySundance import *
from amesosSolver import solverParams


# -------------- set up cell filters ------------------------
class LeftPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt) < 1.0e-10);

class RightPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt-1.0) < 1.0e-10);



# ----------- block triangular solver --------------------

class MyBlockSolver :

  def __init__(self) :
    self.solver_ = buildSolver(solverParams)
    
  def solve(self, A, b, x) :

    A00 = A.getBlock(0,0);
    A11 = A.getBlock(1,1);
    A01 = A.getBlock(0,1);

    b0 = b.getBlock(0)
    b1 = b.getBlock(1)

    x = b.copy()
    x0 = x.getBlock(0)
    x1 = x.getBlock(1)

    state = self.solver_.solve(A11, b1, x1)

    c = b0 - A01*x1
    state = self.solver_.solve(A00, c, x0)

    x.setBlock(0, x0)
    x.setBlock(1, x1)

    return x


# ------------------------------------------------------------


def main():
  """Poisson example code"""
  vecType = EpetraVectorType()
  nProc = getNProc()
  nx = 10
  mesher  = PartitionedLineMesher(0.0, 1.0, nx*nProc);
  mesh = mesher.getMesh();

  u = UnknownFunction(Lagrange(5), "u");
  v = UnknownFunction(Lagrange(3), "v");
  du = TestFunction(Lagrange(5), "du");
  dv = TestFunction(Lagrange(3), "dv");
  
  x = CoordExpr(0);
  dx = Derivative(0);

  quad = GaussianQuadrature(8)
  
  
  interior = MaximalCellFilter()
  points = DimensionalCellFilter(0)
  leftPt = points.subset(LeftPointPredicate())
  rightPt = points.subset(RightPointPredicate())
  
  eqn = Integral(interior, \
                 (dx*du)*(dx*u) + du*v + (dx*dv)*(dx*v) + x*dv, \
                          quad)
  bc = EssentialBC(leftPt, du*u + dv*v, quad) \
       + EssentialBC(rightPt, du*u + dv*v, quad)

  varBlock = BlockList(Block(du, vecType), Block(dv, vecType))
  unkBlock = BlockList(Block(u, vecType),  Block(v, vecType))

  print 'setting up prob'
  prob = LinearProblem(mesh, eqn, bc, varBlock, unkBlock)

  print 'setting up solver'
  solver = MyBlockSolver()

  soln = prob.solve(solver)
  

  x2 = x*x;
  x3 = x*x2;

  uExact = (1.0/120.0)*x2*x3 - 1.0/36.0 * x3 + 7.0/360.0 * x;
  vExact = 1.0/6.0 * x * (x2 - 1.0);

  vDiff = (soln[1] - vExact)**2.0
  uDiff = (soln[0] - uExact)**2.0
  vDiffDeriv = (dx*(soln[1] - vExact))**2.0
  uDiffDeriv = (dx*(soln[0] - uExact))**2.0

  uError = math.sqrt(uDiff.integral(interior, mesh, quad))
  uDerivError = math.sqrt(uDiffDeriv.integral(interior, mesh, quad))

  vError = math.sqrt(vDiff.integral(interior, mesh, quad))
  vDerivError = math.sqrt(vDiffDeriv.integral(interior, mesh, quad))
  
  print "u error = " , uError
  print "u deriv error = " , uDerivError
  
  print "v error = " , vError
  print "v deriv error = " , vDerivError

  uError = max(uError, uDerivError)
  vError = max(vError, vDerivError)
  error = max(uError, vError)

  print "max err = ", error

  writer = MatlabWriter("BlockTest.dat")
  writer.addMesh(mesh)
  writer.addField("u", soln[1])
  writer.addField("v", soln[0])
  writer.write()

  tol = 1.0e-11
  passFailTest(error, tol)
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
