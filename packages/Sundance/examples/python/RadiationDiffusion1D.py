#! /usr/bin/env python

import PySundance
import math
from PySundance import *
from noxSolver import solverParams


# Class PointPredicate tests whether a cell is at a specified point
class PointPredicate :
  # constructor accepts a point position against which cells are tested
  def __init__(self, pointPos) :
    self.pointPos = pointPos
  # return true if cell position is at the specified position pointPos
  def evalOp(self, pt) :
    return (math.fabs(pt-self.pointPos) < 1.0e-10);


# Parameters for linear solver
aztecSolverParams = {"Type" : "Aztec",
                     "Method" : "GMRES",
                     "Max Iterations" : 1000,
                     "Tolerance" : 1.0e-10,
                     "Precond" : "Domain Decomposition",
                     "Subdomain Solver" : "ILU",
                     "Graph Fill" : 1,
                     "Verbosity" : 4
                     }

# Parameters for nonlinear solver
noxSolverParams = {"NOX Solver" :
                   {"Nonlinear Solver" : "Line Search Based",
                    "Line Search" : {"Method" : "More'-Thuente"},
                    "StatusTest"  : {"Max Iterations" : 20,
                                     "Tolerance" : 1.0e-8},
                    "Linear Solver" : aztecSolverParams
                    }
                   }



# --------------------------------------------------------------------------
#
# Solves the radiation diffusion equation in one dimension
#
# --------------------------------------------------------------------------

def main():

  # Create a mesh on a line
  nx = 40
  mesher  = PartitionedLineMesher(0.0, 1.0, nx*getNProc());
  mesh = mesher.getMesh();

  # Create test and unknown functions, expanded in cubic Lagrange basis
  basis = Lagrange(3)
  u = UnknownFunction(basis);
  v = TestFunction(basis);

  # Create differential operator and coordinate function
  x = CoordExpr(0);
  dx = Derivative(0);

  # Set up a discrete function for the initial guess. 
  vecType = EpetraVectorType()
  discSpace = DiscreteSpace(mesh, basis, vecType)
  projector = L2Projector(discSpace, 1.0 + x)
  u0 = projector.project()

  # Define cell filters for the interior and endpoints. The endpoint filters
  # use the PointPredicate defined above
  interior = MaximalCellFilter()
  points = DimensionalCellFilter(0)
  leftPt = points.subset(PointPredicate(0.0))
  rightPt = points.subset(PointPredicate(1.0))

  # Set up the equation and boundary conditions
  quad = GaussianQuadrature(8)
  bc = EssentialBC(leftPt, v*(u-1.0), quad) \
       + EssentialBC(rightPt, v*(u-2.0), quad)
  eqn = Integral(interior, u*u*u*(dx*v)*(dx*u), quad)

  # Create nonlinear problem
  prob = NonlinearProblem(mesh, eqn, bc, v, u, u0, vecType)

  # Create nonlinear solver
  solver = NOXSolver(noxSolverParams, prob)

  # Solve the problem! The expression u0 will contain the solution
  solver.solve()

  # Check against exact solution
  exactSoln = pow(15.0*x + 1.0, 0.25)
  diff = (u0 - exactSoln)**2

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
