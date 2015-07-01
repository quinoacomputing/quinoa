/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */


#include "Sundance.hpp"
//#include "CopyDiscreteFunction.hpp"

/** 
 * Solves the Navier-Stokes equation in 2D using Taylor-Hood-Elements:
 *
 * du/dt - nu* Delta u + (u*grad) u + grad p = f     in Omega x I
 *                                     div u = 0     in Omega x I
 *  Driven Cavity scenario
 */

// Will be used to select the left, bottom, right, and top boundary
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-2.5) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-0.41) < 1.0e-10;})

REFINE_MESH_ESTIMATE(MeshRefEst , { \
    // this is the refinement function
        if ((cellPos[0] < 0.7) && (cellPos[0] > 0.01) && (cellPos[1] < 0.33) && (cellPos[1] > 0.1) && (cellLevel < 2)) \
                                      return 1;\
                                 else return 0; } ,
  // this is the coarse load estimator
  { if (cellPos[0] < 0.7) {return 5;} else {return 1;} } )

MESH_DOMAIN( MeshDomain , {return true;})


int main(int argc, char** argv)
{
  
  try
	{
      Sundance::init(&argc, &argv);
      int myrank = MPIComm::world().getRank();
 
      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();
      
      /* refinement criterion object */
      RefinementClass refCl1 = new MeshRefEst();
      /* Create one circle */
      ParametrizedCurve circle = new Circle( 0.2 , 0.2 , 0.05 , 1 , 0.00001 );
      ParametrizedCurve box = new Box2D( 0.15 , 0.15 , 0.10 , 0.10 , 1 , 0.00001 );
      /* mesh domain estimation */
      MeshDomainDef meshDom_circl = new CurveDomain( box , Outside_Curve );
      /* curve predicate with the same curve! */
      CellPredicate curveIN = new CellCurvePredicate( box , Inside_Curve );
      /* mesh domain, dummy class */
      MeshDomainDef meshDom = new MeshDomain();

      MeshType meshType = new HNMeshType2D();
      MeshSource mesher = new HNMesher2D(0.0, 0.0, 2.5 , 0.41 , 20 , 6 , meshType , refCl1 , meshDom_circl );
      
      Mesh mesh = mesher.getMesh();

      cout << "Nr Points  "<<mesh.numCells(0) << std::endl;
      cout << "My Rank is :" << myrank << std::endl;

      // Create cell filters
      CellFilter interior = new MaximalCellFilter();
      CellFilter nodes = new DimensionalCellFilter(0);
      CellFilter edges = new DimensionalCellFilter(1);

      // Find left, right, top, bottom and cylinder boundary edges
      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());
      CellFilter curveFilt = edges.subset(curveIN);

       // Create unknown and test functions, discretized using Taylor Hood elements 
      BasisFamily L1 = new Lagrange(1);
      BasisFamily L2 = new Lagrange(2);
      Expr ux = new UnknownFunction(L2, "u_x");
      Expr vx = new TestFunction(L2, "v_x");
      Expr uy = new UnknownFunction(L2, "u_y");
      Expr vy = new TestFunction(L2, "v_y");
      Expr p = new UnknownFunction(L1, "p");
      Expr q = new TestFunction(L1, "q");
      Expr u = List(ux, uy);
      Expr v = List(vx, vy);

      // Create differential operator and coordinate functions 
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      // We need a quadrature rule for doing the integrations 
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);
      
      DiscreteSpace VelSpace(mesh, List(L2,L2), vecType);
      DiscreteSpace VelPreSpace(mesh, List(L2,L2,L1), vecType);

      // iterate from last timestep u_old 
      Expr uo = new DiscreteFunction(VelSpace, 0.0, "uo");
      // actual iterate (u,p) 
      Expr up = new DiscreteFunction(VelPreSpace, 0.0, "up");

      // zero right hand side 
      Expr fx = 0.0;
      Expr fy = 0.0;
      // viscosity parameter 
      Expr nu = new Sundance::Parameter(1.0/50.0);
      Expr h = new CellDiameterExpr();

      // initial value 
      L2Projector projectorupar(VelSpace, List(0.0,0.0)); 
    
      L2Projector projectoruparp(VelPreSpace, List(0.0,0.0,0.0)); 
     
      // Define the Dirichlet BC 
      Expr uInflow = 0.5*((0.41-y)*(y));
      CellFilter walls = bottom + top + curveFilt;
      Expr bc = EssentialBC(walls, v*u, quad4 )
              + EssentialBC(left, vx*(ux-uInflow)+vy*uy, quad4 );

      // Set up the variational form of the Navier Stokes equations 
      // or implicitly is defined the zero Neuman boundary conditon
//                             div u = 0     in Omega x I
      Expr eqnp = Integral(interior,(dx*ux+dy*uy)*q,
                           quad2 );
// nu* Delta u + (u*grad) u + grad p = f     in Omega x I
      Expr eqn = Integral(interior,
		nu*(grad*vx)*(grad*ux) + nu*(grad*vy)*(grad*uy) 
		+ vx*(u*grad)*ux + vy*(u*grad)*uy 
		- p*(dx*vx+dy*vy) - vx*fx - vy*fy
		, quad4); 
	//+ Integral(left + right, q*(-1.0), quad4 );
 

      // Create a Playa NonlinearOperator object 
      NonlinearProblem prob(mesh, eqn+eqnp, bc, List(v, q),
                               List(u, p), up, vecType );

      ParameterXMLFileReader reader("nox-newton.xml");
      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams);

      // Solve the nonlinear system
      SolverState<double> status = prob.solve(solver);

      // Write the field in VTK format
      FieldWriter w = new VTKWriter("NavStok_Chanel_box_stat");
      w.addMesh(mesh);

      // this is what should be passed for vector field plotting
      Expr expr_vector(List(up[0],up[1]));

      w.addField("ux", new ExprFieldWrapper(up[0]));
      w.addField("uy", new ExprFieldWrapper(up[1]));

      // add a vector field , try pnly one at once
      w.addField("vel", new ExprFieldWrapper(expr_vector));

      w.addField("p", new ExprFieldWrapper(up[2]));
      w.write();

      double tol = 1.0e-6;
      Expr errExpr1 = Integral(right , (up[0]-uInflow)*(up[0]-uInflow) , quad4);
      FunctionalEvaluator errInt1(mesh, errExpr1);
      double errorSq = errInt1.evaluate();
      double err_i = std::sqrt(errorSq);
      Out::os() << " error ux =" << err_i << std::endl;
      Sundance::passFailTest( err_i , tol);
      
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
  return(0);
}
