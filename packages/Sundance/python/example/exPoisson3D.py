#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance


import math
from PySundance import *

def main():

    vecType = EpetraVectorType()
    mesher  = ExodusNetCDFMeshReader("cube-coarse.ncdf");
    mesh = mesher.getMesh();
    basis = Lagrange(2)

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
    side1 = faces.labeledSubset(1);
    side2 = faces.labeledSubset(2);
    side3 = faces.labeledSubset(3);
    side4 = faces.labeledSubset(4);
    side5 = faces.labeledSubset(5);
    side6 = faces.labeledSubset(6);

    exactSoln = (x + 1.0)*x - 1.0/4.0;

    eqn = Integral(interior, (grad*v)*(grad*u) +2.0*v, quad2)
    bc = EssentialBC(side4, v*(u-exactSoln), quad4) \
        + EssentialBC(side6, v*(u-exactSoln), quad4);

    prob = LinearProblem(mesh, eqn, bc, v, u, vecType)

    solver = readSolver(searchForFile("aztec.xml"));

    soln = prob.solve(solver)

    discSpace = DiscreteSpace(mesh, basis, vecType)
    proj1 = L2Projector(discSpace, exactSoln)
    proj2 = L2Projector(discSpace, soln - exactSoln)
    proj3 = L2Projector(discSpace, pow(soln - exactSoln, 2.0))
    exactDisc = proj1.project()
    error = proj2.project()
    errorSq = proj3.project()

    w = VTKWriter("Poisson3D")
    w.addMesh(mesh);
    w.addField("soln", soln)
    w.addField("exact soln", exactDisc)
    w.addField("error", error)
    w.addField("errorSq", errorSq)
    w.write();

    diff = (soln - exactSoln)**2.0
    diffDeriv = (grad*(soln - exactSoln))**2.0



    error = math.sqrt(diff.integral(interior, mesh, quad4))
    derivError = math.sqrt(diffDeriv.integral(interior, mesh, quad4))
    print "error = " , error
    print "deriv error = " , derivError

    error = max(error, derivError)

    print "max err = ", error

    tol = 1.0e-12
    passFailTest(error, tol)





# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
    main()
