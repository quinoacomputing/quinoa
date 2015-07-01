#! /usr/bin/env python

# @HEADER

# @HEADER

import math
import setpath
import PySundance


from PySundance import *

class LeftPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt) < 1.0e-10);

class RightPointPredicate :
  def evalOp(self, pt) :
    return (math.fabs(pt-math.pi) < 1.0e-10);

def main():

  vecType = EpetraVectorType()
  nx = 100
  nProc = getNProc()
  mesher  = PartitionedLineMesher(0.0, math.pi, nx*nProc);
  mesh = mesher.getMesh();
  basis = Lagrange(4)

  ur = UnknownFunction(basis, "Re(u)");
  vr = TestFunction(basis, "Re(v)");
  ui = UnknownFunction(basis, "Im(u)");
  vi = TestFunction(basis, "Im(v)");

  I = Complex(Expr(0.0), Expr(1.0));

  u = ur + I*ui
  v = conj(vr + I*vi)
  
  x = CoordExpr(0);
  dx = Derivative(0);
  
  quad = GaussianQuadrature(2)
  quad8 = GaussianQuadrature(8)
  
  interior = MaximalCellFilter()
  points = DimensionalCellFilter(0)
  leftPt = points.subset(LeftPointPredicate())
  rightPt = points.subset(RightPointPredicate())

  print 'forming bcs'
  bc = EssentialBC(leftPt, v*(u-1.0), quad) \
       + EssentialBC(rightPt, v*(u+cosh(x)), quad)
  
  eqn = Integral(interior, (dx*v)*(dx*u) - 2.0*I*v*u, quad)
  
  prob = LinearProblem(mesh, eqn, bc, List(vr, vi), List(ur, ui), vecType)

  solver = readSolver(searchForFile("aztec.xml"))

  uSolnRI = prob.solve(solver)

  uSoln = uSolnRI[0] + I*uSolnRI[1]

  print 'solved problem, checking solution'

  exactSoln = cos(x)*cosh(x) - I*sin(x)*sinh(x)

  print 'exact solution = ', exactSoln

  diff = (uSoln-exactSoln).conj() * (uSoln - exactSoln)
  errSq = diff.integral(interior, mesh, quad8)

  print "errorSq = " , errSq
  error = math.sqrt(errSq)
  print "error = " , error

  tol = 1.0e-8
  passFailTest(error, tol)

  mpiSynchronize()
  print "proc ", getRank(), " is done"
  mpiSynchronize()

  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
