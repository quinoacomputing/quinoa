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
    mesher  = ExodusNetCDFMeshReader("widget1.ncdf")
    mesh = mesher.getMesh()

    # Define cell filters that will identify the maximal cells, the
    # inflow surface, and the outflow surface
    interior = MaximalCellFilter()
    inflow = boundary.labeledSubset(1)
    outflow = boundary.labeledSubset(2)
    hole = boundary.labeledSubset(3)
    wall = boundary.labeledSubset(4)

    # Define the unknown and test functions 
    vx = TestFunction(Lagrange(1), "vx");
    vy = TestFunction(Lagrange(1), "vy");
    ux = UnknownFunction(Lagrange(1), "ux");
    uy = UnknownFunction(Lagrange(1), "uy");
                         

    # Define the quadrature rule to be used
    quad2 = GaussianQuadrature(2)

    # Create differential operators
    dx = Derivative(0);
    dy = Derivative(1);
    dz = Derivative(2);
    grad = List(dx, dy, dz);

    E = 1.0
    nu = 0.25
    lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu)
    mu = E/2.0/(1.0 + nu)

    # Create some expressions that will be used in the outflow BC
    x = CoordExpr(0)
    x = CoordExpr(1)
    
    # Set up the weak equation
    strain = List(dx*ux, dy*uy, dx*uy + dy*ux);
    varStrain = List(dx*vx, dy*vy, dx*vy + dy*vx);

    D = List(List(lambda+2.0*mu,    lambda,             0.0), 
             List(lambda,           lambda+2.0*mu,      0.0),
             List(0.0,              0.0,                mu));
    eqn = Integral(interior, varStrain*(D*strain), quad2)

    # Set up the Dirichlet BC
    eps = 0.05
    bc = EssentialBC(outflow, ux*vx + uy*vy, quad2) + \
         EssentialBC(inflow, ux*vx + uy*vy, quad2) + \
         EssentialBC(wall, ux*vx + uy*vy, quad2) + \
         EssentialBC(hole, vx*(ux - eps*(x-0.6) + eps*(y-0.6)*x
                               + vy*(uy - eps*(x-0.6)/2.0), quad2)
         
    


    # Put it all together to create the linear problem object
    prob = LinearProblem(mesh, eqn, bc, List(vx,vy), List(ux,uy), vecType)


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

    # Finally, we can write the velocity and velocity potential to a VTK file
    writer = VTKWriter("MovedMeshD")
    writer.addMesh(mesh)
    writer.addField("u_x", soln[0])
    writer.addField("u_y", soln[1])

    # all done!

    

if __name__ == "__main__":
    main()

    







    
    
