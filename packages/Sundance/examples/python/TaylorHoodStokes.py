#! /usr/bin/env python

# @HEADER

# @HEADER


import PySundance


import math
from PySundance import *

def main():
    vecType = EpetraVectorType()

    # read the mesh from an NCDF file
    mesher  = ExodusNetCDFMeshReader("square-0.01.ncdf");
    mesh = mesher.getMesh();

    # Define the unknown and test functions for velocity and pressure.
    # We will use linear interpolation for all variables, which will
    # require the addition of a stabilization term.
    pOrder = 2
    pressure_basis = FIATLagrange(pOrder)
    velocity_basis = FIATLagrange(pOrder+1)
    ux = UnknownFunction(velocity_basis,'ux');
    vx = TestFunction(velocity_basis,'vx');
    uy = UnknownFunction(velocity_basis,'uy');
    vy = TestFunction(velocity_basis,'vy');
    p = UnknownFunction(pressure_basis,'p');
    q = TestFunction(pressure_basis,'q');

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

    nodes = DimensionalCellFilter(0)
    pinpoint = nodes.labeledSubset(5)
    
    quad = GaussianQuadrature(2*pOrder)
    quadHigh = GaussianQuadrature(2*(pOrder+1))

    eqn = Integral( interior ,
                    (grad * vx) * (grad * ux) \
                    + (grad * vy) * (grad * uy) \
                    - p * (dx * vx + dy * vy ) \
                    + q * (dx * ux + dy * uy ) , quad )

    bc = EssentialBC(left, vx*ux + vy*uy, quadHigh) \
         + EssentialBC(right, vx*ux + vy*uy, quadHigh) \
         + EssentialBC(top, vx*(ux-1.0) + vy*uy, quadHigh) \
         + EssentialBC(bottom, vx*ux + vy*uy, quadHigh) \
         + EssentialBC(pinpoint, q*p, quad);


#    vecBasis = BasisList(velocity_basis, velocity_basis, pressure_basis);
#    discSpace = DiscreteSpace(mesh, vecBasis, vecType);
#    u0 = DiscreteFunction(discSpace, 0.0);

    linSolver = readSolver("aztec.xml");

    print 'prob ctor'
    prob = LinearProblem(mesh, eqn, bc, List(vx,vy,q), List(ux,uy,p), vecType)
    #solver = NOXSolver(solverParams, prob)

    print 'solving...'
    u0 = prob.solve(linSolver)
    print '...done solve'

    # set up streamfunction calculation
    basis = FIATLagrange(1)
    psi = UnknownFunction(basis, "psi")
    varPsi = TestFunction(basis, "varPsi")

    bdry = BoundaryCellFilter()

    ux0 = u0[0]
    uy0 = u0[1]

    curlU0 = dx*uy0 - dy*ux0
    psiEqn = Integral(interior, (grad*psi)*(grad*varPsi) + varPsi*curlU0, quad)
    psiBC = EssentialBC(bdry, varPsi*psi, quad)
    psiProb = LinearProblem(mesh, psiEqn, psiBC, varPsi, psi, vecType)
    


    psi0 = psiProb.solve(linSolver)

    w = VTKWriter("Stokes")
    w.addMesh(mesh);
    w.addField("psi", psi0)
    w.write();




# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
    main()
import PySundance
import math

