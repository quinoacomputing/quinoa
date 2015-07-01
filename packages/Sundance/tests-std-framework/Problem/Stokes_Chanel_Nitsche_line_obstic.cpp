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

/** 
 * Solves the Stokes equation in 2D using Taylor-Hood-Elements:
 *
 * Delta u + grad p = f     in Omega x I
 *            div u = 0     in Omega x I
 * Using the Nitsche method to impose the boundary
 */

// Will be used to select the left, bottom, right, and top boundary
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-2.2) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-0.41) < 1.0e-10;})

CELL_PREDICATE(RightBottomPointTest, { return ((fabs(x[1]) < 1.0e-10)&&(fabs(x[0]-2.2) < 1.0e-10)); })

// define here the refinement level
#define REFINE_LEVEL 1

REFINE_MESH_ESTIMATE(MeshRefEst , { \
    // this is the refinement function , if the returned value > 0 then makes a refinement in the given cellLevel , 
    // the function is evaluated per cell once
        if ((cellPos[0] > 2.0) && (cellLevel < REFINE_LEVEL )) \
                                      return 1;\
                                 else return 0; } ,
  // this is the coarse load estimator , only for the parallelization of the mesh
  { if (cellPos[0] < 0.7) {return 5;} else {return 1;} } )

// dummy mesh domain estimator
MESH_DOMAIN( MeshDomain , {return true;})


int main(int argc, char** argv)
{
  
  try
	{
      Sundance::init(&argc, &argv);
      //int np = MPIComm::world().getNProc();
      int myrank = MPIComm::world().getRank();
 
      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();
      
      /* refinement criterion object */
      RefinementClass refCl1 = new MeshRefEst();
      /* Create one circle */
      double offsPoint = 0.153;
      ParametrizedCurve curve = new Line2D( -1.0 ,  2.2+0.41-offsPoint , 1 , 1e-8 );
            
      /* Obiject needed for curve integral calculation */
      ParametrizedCurve curveIntegral = new ParamCurveIntegral(curve);
      
      /* mesh domain estimation , cells which are inside the curve will be deleted*/
      MeshDomainDef meshDom_circl = new CurveDomain( curve , Outside_Curve );

      /* mesh domain, dummy class , use this if you need all the cells from the mesh (even those in the )*/
      MeshDomainDef meshDom = new MeshDomain();

      // create the mesh
      MeshType meshType = new HNMeshType2D();
      MeshSource mesher = new HNMesher2D(0.0, 0.0, 2.2 , 0.41 , 51 , 11 , meshType , refCl1 , meshDom_circl );
      Mesh mesh = mesher.getMesh();
      ParametrizedCurve curve_line = curve.getPolygon(mesh, 0.01); 
      curve_line.writeToVTK("p_line.vtk");

      cout << "Nr Points  "<<mesh.numCells(0) << endl;
      cout << "My Rank is :" << myrank << endl;

      // Create cell filters
      CellFilter interior = new MaximalCellFilter();
      CellFilter nodes = new DimensionalCellFilter(0);
      CellFilter edges = new DimensionalCellFilter(1);

      // Find left, right, top, bottom and cylinder boundary edges
      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());
      CellFilter right_bottom = nodes.subset(new RightBottomPointTest());
     
      
      /* the cell predicate for different integrations on curve */
      CellPredicate curveIN = new CellCurvePredicate( curve , Inside_Curve );
      CellPredicate curveOUT = new CellCurvePredicate( curve , Outside_Curve );
      CellPredicate curveON = new CellCurvePredicate( curve , On_Curve );
      
      CellFilter InCircle = interior.subset(curveIN);
      CellFilter OutsideCircle = interior.subset(curveOUT);
      CellFilter OnCircle = interior.subset(curveON);
      CellFilter newWall = edges.subset(curveIN);

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
      QuadratureFamily quad2 = new GaussianQuadrature(4);
      QuadratureFamily quad4 = new GaussianQuadrature(4);
      QuadratureFamily quad_c = new GaussianQuadrature(6);
      QuadratureFamily quad_adaptive = new GaussLobattoQuadrature(5);//new TrapesoidQuadrature(100);//new FeketeQuadrature(10); //new GaussLobattoQuadrature(4);
      
      // zero right hand side 
      Expr fx = 0.0;
      Expr fy = 0.0;
      // viscosity parameter 
      Expr nu = new Sundance::Parameter(1.0/50.0);
      Expr h = new CellDiameterExpr();
      // one could use direct values
      Expr nx = new CurveNormExpr(0); //0.707106781186547; //new CurveNormExpr(0);
      Expr ny = new CurveNormExpr(1); //0.707106781186547; //new CurveNormExpr(1);
      double ga1=10;
      double ga2=10;
      // the fixes Dirichlet values on the circle (cylinder)
      double velx = 0.0;
      double vely = 0.0;
     
      // Define the Dirichlet BC 
      Expr uInflow = 0.5*((0.41-y)*(y));
      CellFilter walls = bottom + top;
      Expr bc = EssentialBC(walls, v*u, quad4 )
		+ EssentialBC(left, vy*uy + (ux-uInflow)*vx, quad4 )
		+ EssentialBC( right_bottom , p*q , quad2 );

      // The equation
      Expr eqn =
      // Stokes on the cells which are completly in the computational domain (the 1/Re actually is not necesary)
       Integral( OutsideCircle , nu*(grad*vx)*(grad*ux) +
                              nu*(grad*vy)*(grad*uy) 
                            - p*(dx*vx+dy*vy) + (dx*ux+dy*uy)*q , quad4)
       + Integral(OnCircle, nu*(grad*vx)*(grad*ux) + nu*(grad*vy)*(grad*uy) 
                            - p*(dx*vx+dy*vy) + (dx*ux+dy*uy)*q ,
                  quad_adaptive , curve)

      // Now the Nitsche formula for the boundary condition on the cylinder    
       + Integral(OnCircle, -nu*(dx*ux*nx+dy*ux*ny)*vx-nu*(dx*uy*nx+dy*uy*ny)*vy
                           +p*(nx*vx+ny*vy)
                -nu*(dx*vx*nx+dy*vx*ny)*(ux-velx)-nu*(dx*vy*ny+dy*vy*nx)*(uy-vely)
                           -q*(nx*(ux-velx)+ny*(uy-vely)),
                           quad_c , curveIntegral )
       + Integral(OnCircle, nu*ga1/h*((ux-velx)*vx+(uy-vely)*vy)
                           +ga2/h*((ux-velx)*nx+(uy-vely)*ny)
                                 *(vx*nx+vy*ny)
                           , quad_c , curveIntegral);

      // set up the linear problem 
      LinearProblem prob(mesh, eqn, bc, List(vx, vy , q), 
                         List(ux ,uy ,p), vecType);
  
      ParameterXMLFileReader reader("amesos.xml");
      ParameterList solverParams = reader.getParameters();
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr up = prob.solve(solver);

      /* Calculate the the pressure gradient */
      DiscreteSpace discreteSpace(mesh, List( L1 , L1 ), vecType);
      L2Projector projector(discreteSpace, grad*up[2]);
      Expr p_grad = projector.project();
      
      // Write the field in VTK format
      FieldWriter w = new VTKWriter("Stokes_Chanel_Nitsche_line_obst");
      w.addMesh(mesh);

      // this is what should be passed for vector field plotting
      Expr expr_vector(List(up[0],up[1]));
      Expr expr_p_grad(List(p_grad[0],p_grad[1]));
      w.addField("ux", new ExprFieldWrapper(up[0]));
      w.addField("uy", new ExprFieldWrapper(up[1]));
      // add a vector field to the plot
      w.addField("vel", new ExprFieldWrapper(expr_vector));
      w.addField("p", new ExprFieldWrapper(up[2]));
      w.addField("p_grad", new ExprFieldWrapper(expr_p_grad));
      w.addField("p_grad_x", new ExprFieldWrapper(p_grad[0]));
      w.addField("p_grad_y", new ExprFieldWrapper(p_grad[1]));
      w.write();
    
      Expr curveExp = Integral(OnCircle, 1 ,quad_c , curveIntegral );
      FunctionalEvaluator curveInt( mesh , curveExp);
      double RcurveInt = curveInt.evaluate();
      Out::os() << "curveInt int 1 dc = " << RcurveInt << std::endl;
      Sundance::passFailTest( (RcurveInt - std::sqrt(2*offsPoint*offsPoint)) , 1e-6);
      
      Expr curveExpX = Integral(OnCircle, x ,quad_c , curveIntegral );
      FunctionalEvaluator curveIntX( mesh , curveExpX);
      double RcurveIntX = curveIntX.evaluate();
      Out::os() << "curveIntX = " << RcurveIntX << std::endl;
      double lx = offsPoint;
      double analyticSolCX = (std::sqrt(2.0)/2.0)*(2.2*2.2 - (2.2-lx)*(2.2-lx));
      Out::os() << "curveIntX(analytic) int x dc = " << analyticSolCX << std::endl;
      Sundance::passFailTest( (RcurveIntX-analyticSolCX) , 1e-6);
      
      Expr curveExpY = Integral(OnCircle, y ,quad_c , curveIntegral );
      FunctionalEvaluator curveIntY( mesh , curveExpY);
      double RcurveIntY = curveIntY.evaluate();
      Out::os() << "curveIntY = " << RcurveIntY << std::endl;
      double analyticSolCY = (std::sqrt(2.0)/2.0)*(0.41*0.41 - (0.41-lx)*(0.41-lx));
      Out::os() << "curveIntY(analytic) int y dc = " << analyticSolCY << std::endl;
      Sundance::passFailTest( (RcurveIntY-analyticSolCY) , 1e-6);
      
      Expr areaExp = Integral(OutsideCircle , 1 , quad2 ) + Integral(OnCircle , 1 , quad_adaptive , curve ); // watchMe
      FunctionalEvaluator areInt( mesh , areaExp);
      double area = areInt.evaluate();
      Out::os() << "Area = " << area << std::endl;
      Sundance::passFailTest( (area-(0.41*2.2-0.5*offsPoint*offsPoint)) , 1e-6);
      
      Expr funcExp1 = Integral(OutsideCircle , x*y , quad4 );
      Expr funcExp2 = Integral(OnCircle , x*y , quad_adaptive , curve ); // watchMe
      FunctionalEvaluator funcInt1( mesh , funcExp1);
      FunctionalEvaluator funcInt2( mesh , funcExp2);
      double funcI1 = funcInt1.evaluate();
      double funcI2 = funcInt2.evaluate();
      Out::os() << " int_{whole cell} (x*y) = " << funcI1 << std::endl;
      Out::os() << " int_{cutt cell} (x*y) = " << funcI2 << std::endl;
      Out::os() << " int_{omega} (x*y) = " << funcI1+funcI2 << std::endl;
      // here we calculate the analytical solution of the integral
      double b = 2.2+0.41-lx;
      double i1 = ( ((0.41*0.41) / 2.0) * (((2.2-lx)*(2.2-lx))/2.0));
      double x2 = 2.2;
      double F2 = (0.5*( (1.0/2.0)*b*b*x2*x2 - (2.0/3)*b*x2*x2*x2 + x2*x2*x2*x2/4.0));
      double x1 = 2.2-lx;
      double F1 = (0.5*( (1.0/2.0)*b*b*x1*x1 - (2.0/3.0)*b*x1*x1*x1 + x1*x1*x1*x1/4.0));
      double totalInt = i1 + (F2-F1);
      Out::os() << "(Analytic) int_{omega} (x*y)  = " << totalInt << std::endl;
      Sundance::passFailTest( (funcI1+funcI2-totalInt) , 1e-6); //0.196237948515625 -> the are with 0.135
      
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
  return(0);
}
