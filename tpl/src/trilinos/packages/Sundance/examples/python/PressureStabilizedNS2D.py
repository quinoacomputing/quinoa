#! /usr/bin/env python

# @HEADER

# @HEADER

import setpath
import PySundance


import math
from PySundance import *

from noxSolver import solverParams

def main():
    

    vecType = EpetraVectorType()

    # read the mesh from an NCDF file
    mesher  = ExodusNetCDFMeshReader("../../../examples-tutorial/meshes/square-0.01.ncdf");
    mesh = mesher.getMesh();

    # Define the unknown and test functions for velocity and pressure.
    # We will use linear interpolation for all variables, which will
    # require the addition of a stabilization term.
    basis = Lagrange(1)
    ux = UnknownFunction(basis,'ux');
    vx = TestFunction(basis,'vx');
    uy = UnknownFunction(basis,'uy');
    vy = TestFunction(basis,'vy');
    p = UnknownFunction(basis,'p');
    q = TestFunction(basis,'q');

    # Group the velocities into a vector-valued expression
    u = List(ux, uy)
    v = List(vx, vy)

    # Create differential operators
    dx = Derivative(0);
    dy = Derivative(1);
    grad = List(dx, dy);

    interior = MaximalCellFilter()
    edges = DimensionalCellFilter(1)
    bottom = edges.labeledSubset(1);
    right = edges.labeledSubset(2);
    top = edges.labeledSubset(3);
    left = edges.labeledSubset(4);
    
    h = 0.02;
    beta = 0.08;
    reynolds = Parameter(1.0);

    quad2 = GaussianQuadrature(2)
    quad4 = GaussianQuadrature(4)


    eqn = Integral(interior, (grad*vx)*(grad*ux)  \
                   + (grad*vy)*(grad*uy)  - p*(dx*vx+dy*vy) \
                   + beta*h*h*(grad*q)*(grad*p) + q*(dx*ux+dy*uy), \
                   quad2) \
                   + Integral(interior, reynolds*(vx*(u*grad)*ux) \
                              + reynolds*(vy*(u*grad)*uy), quad2);

    bc = EssentialBC(left, vx*ux + vy*uy, quad2) \
         + EssentialBC(right, vx*ux + vy*uy, quad2) \
         + EssentialBC(top, vx*(ux-1.0) + vy*uy, quad2) \
         + EssentialBC(bottom, vx*ux + vy*uy, quad2);

    vecBasis = BasisList(basis, basis, basis);
    discSpace = DiscreteSpace(mesh, vecBasis, vecType);
    u0 = DiscreteFunction(discSpace, 0.0);
    
    prob = NonlinearProblem(mesh, eqn, bc, List(vx,vy,q), List(ux,uy,p), u0, vecType)
    solver = NOXSolver(solverParams, prob)

    numReynolds = 20;
    finalReynolds = 500.0;


    # set up streamfunction calculation
    psi = UnknownFunction(basis, "psi")
    varPsi = TestFunction(basis, "varPsi")

    bdry = BoundaryCellFilter()

    ux0 = u0[0]
    uy0 = u0[1]

    curlU0 = dx*uy0 - dy*ux0
    psiEqn = Integral(interior, (grad*psi)*(grad*varPsi) + varPsi*curlU0, quad4)
    psiBC = EssentialBC(bdry, varPsi*psi, quad2)
    psiProb = LinearProblem(mesh, psiEqn, psiBC, varPsi, psi, vecType)
    


    linSolver = readSolver("../../../tests-std-framework/Problem/aztec.xml");
    
    for i in range(numReynolds+1) :
        Re = i*finalReynolds/numReynolds;
        reynolds.setParameterValue(Re);
        print 'doing Reynolds=', Re
        solver.solve()
        psi0 = psiProb.solve(linSolver)
        w = VTKWriter("NS2D-" + str(Re))
        w.addMesh(mesh);
        w.addField("ux", u0[0])
        w.addField("uy", u0[1])
        w.addField("p", u0[2])
        w.addField("psi", psi0)
        w.write();




# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
    main()
