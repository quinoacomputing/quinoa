#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance


import math
from PySundance import *


def runtest(nMesh, order):

  vecType = EpetraVectorType()
  npx = 1
  npy = getNProc()

  mesher  = PartitionedRectangleMesher(0.0, 1.0, nMesh, npx,
                                       0.0, 1.0, nMesh, npy);
  mesh = mesher.getMesh();
  basis = Lagrange(order)

  u = UnknownFunction(basis, "u");
  v = TestFunction(basis, "v");
  x = CoordExpr(0);
  y = CoordExpr(1);
  dx = Derivative(0);
  dy = Derivative(1);
  grad = List(dx, dy);

  quad = GaussianQuadrature(2*order)
  quad6 = GaussianQuadrature(6)

  interior = MaximalCellFilter()
  bdry = BoundaryCellFilter()

  pi = 4.0 * math.atan(1.0)
  eqn = Integral(interior, (grad*v)*(grad*u), quad)\
        + Integral(interior, -v*2.0*pi*pi*order*order*sin(order*pi*x)*sin(order*pi*y), quad)

  bc = EssentialBC(bdry, v*u, quad)

  prob = LinearProblem(mesh, eqn, bc, v, u, vecType)


  solver = readSolver(searchForFile("aztec.xml"));


  soln = prob.solve(solver)

  exactSoln = sin(order*pi*x)*sin(order*pi*y)
  diff = (soln - exactSoln)**2.0

  error = math.sqrt(diff.integral(interior, mesh, quad))


  return error


def linefit(data) :
  sx = 0.0;
  sy = 0.0;
  sxy = 0.0;
  sxx = 0.0;

  for (x,y) in (data) :
    sx += x
    sxx += x*x
    sy += y
    sxy += x*y

  N = len(data)

  D = N*sxx - sx*sx
  a = (sxx*sy - sx*sxy)/D
  b = (N*sxy - sx*sy)/D

  return (a,b)


def main() :

  ok = 1
  
  sizes = (4, 8, 16, 24, 32, 48, 64, 96)

  print sizes

  for order in range(1,2):

    data = []

    for n in sizes:

      print 'order=%d n=%d' % (order,n)
      err = runtest(n, order)
      print 'error = ' , err
      data = data + [(math.log(1.0/n), math.log(err))] 

    (intercept, slope) = linefit(data)
      
    print "order=", order, " slope = " , slope
      
    if math.fabs(slope-(order+1)) > 0.2 :
      ok = 0
      print "FAILED"
    else:
      print "PASSED"
      
  if (ok != 1) :
    print "FAILED: failures detected"
  else:
    print "all orders PASSED"


# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
