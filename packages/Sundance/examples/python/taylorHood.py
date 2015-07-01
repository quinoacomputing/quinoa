#!/usr/bin/env python



import setpath
import PySundance
import math
from PySundance import *


from noxSolver import solverParams

def main() : 

    vecType = EpetraVectorType()

    # read the mesh from an NCDF file
    mesher  = ExodusNetCDFMeshReader("../../examples-tutorial/meshes/square-0.01.ncdf");
    mesh = mesher.getMesh();

    # Define the unknown and test functions for velocity and pressure.
    # We will use linear interpolation for all variables, which will
    # require the addition of a stabilization term.
    Pk = Lagrange(1)
    Pk1 = Lagrange(2)
    ux = UnknownFunction(Pk1,'ux');
    vx = TestFunction(Pk1,'vx');
    uy = UnknownFunction(Pk1,'uy');
    vy = TestFunction(Pk1,'vy');
    p = UnknownFunction(Pk,'p');
    q = TestFunction(Pk,'q');

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

    
    reynolds = Parameter(1.0);

    quad2 = GaussianQuadrature(2)
    quad4 = GaussianQuadrature(4)

    eqn = Integral(interior, (grad*vx)*(grad*ux)  \
                   + (grad*vy)*(grad*uy)  - p*(dx*vx+dy*vy) + q*(dx*ux+dy*uy), quad4) \
                   + Integral(interior, reynolds*(vx*(u*grad)*ux) \
                              + reynolds*(vy*(u*grad)*uy), quad4);

    bc = EssentialBC(left, vx*ux + vy*uy, quad2) \
         + EssentialBC(right, vx*ux + vy*uy, quad2) \
         + EssentialBC(top, vx*(ux-1.0) + vy*uy, quad2) \
         + EssentialBC(bottom, vx*ux + vy*uy, quad2);

    vecBasis = BasisList(Pk1, Pk1, Pk);
    discSpace = DiscreteSpace(mesh, vecBasis, vecType);
    u0 = DiscreteFunction(discSpace, 0.0);

    print 'prob ctor'
    prob = NonlinearProblem(mesh, eqn, bc, List(vx,vy,q), List(ux,uy,p), u0, vecType)

    print 'solver ctor'
    solver = NOXSolver(solverParams, prob)

    numReynolds = 10;
    finalReynolds = 200.0;

    print 'starting continuation loop'
    
    for i in range(numReynolds+1) :
        Re = i*finalReynolds/numReynolds;
        reynolds.setParameterValue(Re);
        print 'doing Reynolds=', Re
        solver.solve()
        w = VTKWriter("NS2D-" + str(Re))
        w.addMesh(mesh);
        w.addField("ux", u0[0])
        w.addField("uy", u0[1])
        w.addField("p", u0[2])
        w.write();



    

    
if __name__ == "__main__":
    main()
