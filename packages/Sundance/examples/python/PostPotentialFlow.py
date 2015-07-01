#! /usr/bin/env python


import PySundance
import math
from PySundance import *


#############################################################################
#
# Potential flow example
#
#############################################################################

def main():
    vecType = EpetraVectorType()

    # Read the mesh from an Exodus NCDF file
    mesher  = ExodusNetCDFMeshReader("../meshes/post.ncdf")
    mesh = mesher.getMesh()

    # Define cell filters that will identify the maximal cells, the
    # inflow surface, and the outflow surface
    interior = MaximalCellFilter()
    boundary = BoundaryCellFilter()
    inflow = boundary.labeledSubset(1)
    outflow = boundary.labeledSubset(2)

    # Define the unknown and test functions 
    phi = UnknownFunction(Lagrange(1), "phi")
    phiHat = TestFunction(Lagrange(1), "phiHat")
                         

    # Define the quadrature rule to be used
    quad2 = GaussianQuadrature(2)

    # Create differential operators
    dx = Derivative(0);
    dy = Derivative(1);
    dz = Derivative(2);
    grad = List(dx, dy, dz);

    # Create some expressions that will be used in the outflow BC
    x = CoordExpr(0)
    L = 1.0;

    # Set up the weak equation
    eqn = Integral(interior, (grad*phiHat)*(grad*phi), quad2) \
          + Integral(inflow, phiHat*(x-phi)/L, quad2)

    # Set up the Dirichlet BC
    bc = EssentialBC(outflow, phiHat*phi, quad2)


    # Put it all together to create the linear problem object
    prob = LinearProblem(mesh, eqn, bc, phiHat, phi, vecType)


    # set up a linear solver.
    solverParams = ParameterList({"Linear Solver" : 
                                  {"Type" : "TSF",
                                   "Method" : "BICGSTAB",
                                   "Max Iterations" : 1000,
                                   "Tolerance" : 1.0e-12,
                                   "Precond" : "ILUK",
                                   "Graph Fill" : 1,
                                   "Verbosity" : 4
                                   }
                                  })

    linSolver = buildSolver(solverParams)

    # solve the problem!!!
    soln = prob.solve(linSolver)

    # We've solved the problem; now we want to visualize the result. We'd
    # like to look at the velocity field as well as the velocity potential;
    # however, because the basis for the potential is not C1 the velocity
    # field can't be represented exactly in the same discrete space
    # as the potential. What to do? We'll project the velocity field
    # onto a discrete space of continuous functions so we can visualize it.
    discreteSpace = DiscreteSpace(mesh, BasisList(Lagrange(1),
                                                  Lagrange(1),
                                                  Lagrange(1)), vecType)

    projector = L2Projector(discreteSpace, grad*soln)
    velocity = projector.project()

    # Finally, we can write the velocity and velocity potential to a VTK file
    writer = VTKWriter("Post3D")
    writer.addMesh(mesh)
    writer.addField("phi", soln)
    writer.addField("u_x", velocity[0])
    writer.addField("u_y", velocity[1])
    writer.addField("u_z", velocity[2])

    # all done!

    

if __name__ == "__main__":
    main()

    







    
    
