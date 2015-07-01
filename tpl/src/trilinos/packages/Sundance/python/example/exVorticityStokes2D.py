#!/usr/bin/env python

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
      return (math.fabs(y-1.0) < 1.0e-10);

def main():
  
  skipTimingOutput();
  
  vecType = EpetraVectorType()
  npx = 1
  npy = getNProc()
  ny = 4
  nx = 4
  mesher  = PartitionedRectangleMesher(0.0, 1.0, nx, npx,
                                       0.0, 1.0, ny/npy, npy);
  mesh = mesher.getMesh();

  basis = Lagrange(1)

  psi = UnknownFunction(basis);
  vPsi = TestFunction(basis);
  omega = UnknownFunction(basis);
  vOmega = TestFunction(basis);
  x = CoordExpr(0);
  y = CoordExpr(1);
  dx = Derivative(0);
  dy = Derivative(1);
  grad = List(dx, dy);
  
  quad2 = GaussianQuadrature(2)
  
  interior = MaximalCellFilter()
  edges = DimensionalCellFilter(1)
  left = edges.subset(LeftPointPredicate())
  right = edges.subset(RightPointPredicate())
  top = edges.subset(TopPointPredicate())
  bottom = edges.subset(BottomPointPredicate())

  eqn = Integral(interior, (grad*vPsi)*(grad*psi) \
                 + (grad*vOmega)*(grad*omega) + vPsi*omega, quad2) \
                 + Integral(top, -1.0*vPsi, quad2)
  
  bc = EssentialBC(bottom, vOmega*psi, quad2) \
       + EssentialBC(top, vOmega*psi, quad2) \
       + EssentialBC(left, vOmega*psi, quad2) \
       + EssentialBC(right, vOmega*psi, quad2)
  
  prob = LinearProblem(mesh, eqn, bc, List(vPsi, vOmega), 
                       List(psi, omega), vecType)

  
  solver = readSolver(searchForFile("bicgstab.xml"));

  soln = prob.solve(solver)


  writer = VTKWriter("VorticityStokes2D");
  writer.addMesh(mesh)
  writer.addField("psi", soln[0])
  writer.addField("omega", soln[1])
  writer.write()
  
  
  # As a check, we integrate the vorticity over the domain. By
  # Stokes' theorem this should be equal to the line integral
  # of the velocity around the boundary.
  totalVorticity = soln[1].integral(interior, mesh, quad2)
  error = math.fabs(1.0 - totalVorticity)

  
  
  tol = 1.0e-10
  passFailTest(error, tol)



  
  
# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
