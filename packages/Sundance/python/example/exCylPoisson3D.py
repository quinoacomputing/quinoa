#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance


import math
from PySundance import *

import setpath
from PySundance import *



aztecSolverDict = {"Linear Solver" : 
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

solverParams = ParameterList(aztecSolverDict);


def main():

    vecType = EpetraVectorType()
    mesher  = TriangleMeshReader("ring3D-0.008a");

    print 'reading mesh...'
    mesh = mesher.getMesh();
    basis = Lagrange(1)

    u = UnknownFunction(basis, "u");
    v = TestFunction(basis, "v");
    x = CoordExpr(0);
    y = CoordExpr(1);
    z = CoordExpr(2);
    dx = Derivative(0);
    dy = Derivative(1);
    dz = Derivative(2);
    grad = List(dx, dy, dz);

    quad2 = GaussianQuadrature(2)
    quad4 = GaussianQuadrature(4)

    interior = MaximalCellFilter()
    faces = DimensionalCellFilter(2)
    inner = faces.labeledSubset(1);
    outer = faces.labeledSubset(2);

    rSq = x*x + y*y
    exactSoln = log(rSq)

    eqn = Integral(interior, (grad*v)*(grad*u), quad2)
    bc = EssentialBC(inner, v*(u-exactSoln), quad4) \
        + EssentialBC(outer, v*(u-exactSoln), quad4);

    prob = LinearProblem(mesh, eqn, bc, v, u, vecType)

    print 'solving problem...'
    solver = buildSolver(solverParams)
    soln = prob.solve(solver)

    print 'discretizing exact soln and error norms...'
    discSpace = DiscreteSpace(mesh, basis, vecType)
    proj3 = L2Projector(discSpace, pow(soln - exactSoln, 2.0))
    procField = DiscreteFunction(discSpace, getRank())
    errorSq = proj3.project()

    print 'writing viz for exact soln and error norms...'
    w = VTKWriter("CylPoisson3D")
    w.addMesh(mesh);
    w.addField("soln", soln)
    w.addField("procID", procField)
    w.addField("errorSq", errorSq)
    w.write();

    diff = (soln - exactSoln)**2.0
    diffDeriv = (grad*(soln - exactSoln))**2.0


    error = math.sqrt(diff.integral(interior, mesh, quad4))
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
