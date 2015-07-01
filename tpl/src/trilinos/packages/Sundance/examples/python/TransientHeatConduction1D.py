#!/usr/bin/env python

import PySundance
import math
from PySundance import *

import setpath
from PySundance import *

# Class PointPredicate tests whether a cell is at a specified point
class PointPredicate :
  # constructor accepts a point position against which cells are tested
  def __init__(self, pointPos) :
    self.pointPos = pointPos
  # return true if cell position is at the specified position pointPos
  def evalOp(self, pt) :
    return (math.fabs(pt-self.pointPos) < 1.0e-10);

# define solver parameters:
aztecSolverParams = {"Linear Solver" : 
                     {"Type" : "Aztec",
                      "Method" : "GMRES",
                      "Max Iterations" : 1000,
                      "Tolerance" : 1.0e-12,
                      "Precond" : "Domain Decomposition",
                      "Subdomain Solver" : "ILU",
                      "Graph Fill" : 1,
                      "Verbosity" : 4
                      }
                     }


def main():
    # Specify that we'll do linear algebra with Epetra objects
    vecType = EpetraVectorType()


    # Mesh [0.0, 1.0] with NX elements per processor
    nProc = getNProc()
    nx = 10
    mesher  = PartitionedLineMesher(0.0, 1.0, nx*nProc);
    mesh = mesher.getMesh();

    # Define test and unknown functions
    basis = Lagrange(2)
    u = UnknownFunction(basis, "u");
    v = TestFunction(basis, "v");

    # Define a discrete function used to represent the previous timestep
    discSpace = DiscreteSpace(mesh, basis, vecType)
    pi = 4.0*math.atan(1.0)
    projector = L2Projector(discSpace, sin(pi*x) + 0.3*sin(2.0*pi*x) )

    # Define some symbolic objects used in specifying the equations
    x = CoordExpr(0);
    dx = Derivative(0);

    # Define cell filters used to locate BCs
    interior = MaximalCellFilter()
    points = DimensionalCellFilter(0)
    leftPt = points.subset(PointPredicate(0.0))

    # Define a quadrature rule
    quadOrder = 2
    quad = GaussianQuadrature(quadOrder)
    bc = EssentialBC(leftPt, v*u, quad)
    eqn = Integral(interior, -(dx*v)*(dx*u) - 2*v, quad)
  
    prob = LinearProblem(mesh, eqn, bc, v, u, vecType)

    solver = buildSolver(aztecSolverParams)

    soln = prob.solve(solver)

    exactSoln = x*(x-2.0)

    diff = (soln - exactSoln)**2.0
    error = math.sqrt(diff.integral(interior, mesh, quad))

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
