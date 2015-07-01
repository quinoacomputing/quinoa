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
  
/* We solve in 2D
   - Laplace u = 3.0
   stationary
   With Dirichlet boundary condition:
             u = const , at the boundary
 */
CELL_PREDICATE(LeftPointTest, {return fabs(x[0] - 0.0) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1] - 0.0) < 1.0e-10;})
CELL_PREDICATE(FrontPointTest, {return fabs(x[2]-0.0) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})
CELL_PREDICATE(BackPointTest, {return fabs(x[2]-1.0) < 1.0e-10;})

REFINE_MESH_ESTIMATE(MeshRefEst , { return 0; } , {return 1;} )
MESH_DOMAIN( MeshDomain , {return true;})


REFINE_MESH_ESTIMATE(MeshRefEst1 , { \
	if ( ( (cellPos[0] > 0.5) && (cellPos[0] < 1.0) && (cellPos[1] > 0.5) && (cellPos[1] < 1.0) &&  
	     (cellPos[2] > 0.5) && (cellPos[2] < 1.0) && (cellLevel < 2)) ) 
                  return 1;\
                  else return 0; } , {return 1;} )
int main(int argc, char** argv)
{
  try
  {
      Sundance::init(&argc, &argv);
      //int np = MPIComm::world().getNProc();
      int myrank = MPIComm::world().getRank();

      // We will do our linear algebra using Epetra 
      VectorType<double> vecType = new EpetraVectorType();  
      
      //ParametrizedCurve curve = new Plane3D( 1.0 , 1.0 , -0.4513 , 1 , 0.00001 );
      //ParametrizedCurve curve = new Sphere( 1.0 , 1.0 , 1.0 , 0.4713 , 1 , 0.00001 );
      //ParametrizedCurve curve = new Plane3D( 1e7 , 1e-8 , 5-1e7 , 1 , 0.00001 );

      Array<Point> pts(4);
      Array<int> ptsIndex(6);
      pts[0] = Point( -0.5 , -0.5 , -0.53 );
      pts[1] = Point(  1.5 , -0.5 , -0.53 );
      pts[2] = Point( -0.5 , 1.5 , 1.47 );
      pts[3] = Point(  1.5 , 1.5 , 1.47 );
      ptsIndex[0] = 0; ptsIndex[1] = 1; ptsIndex[2] = 2;
      ptsIndex[3] = 1; ptsIndex[4] = 3; ptsIndex[5] = 2;
      
      ParametrizedCurve curve = new TriangleSurf3D( pts , ptsIndex , 1.0 , 1e-8 );
      //ParametrizedCurve curve = TriangleSurf3D::importGTSSurface( "sphere5.gts" , 1.0 , 1e-8 , Point(0.0,0.0,0.0));
      //ParametrizedCurve curve = TriangleSurf3D::importGTSSurface( "cube.gts" , 1.0 , 1e-8 , Point(0.0,0.0,0.0));
      //ParametrizedCurve curve = TriangleSurf3D::importSTLSurface( "cube.stl" , 1.0 , 1e-8 , Point(0.5,0.5,0.5));
      //curve.flipDomains();
      
      curve.addNewScalarField( "dispX" , 0.0 );
      curve.addNewScalarField( "dispY" , 1.0 );
      curve.addNewScalarField( "dispZ" , 2.0 );
      
      curve.writeToVTK("surf_3D.vtk");
      
      /* refinement criterion object */
      RefinementClass refCl1 = new MeshRefEst1();
      RefinementClass refCl = new MeshRefEst();
      RefinementClass refGeom = new GeometryRefinement(curve,2);

      /* mesh domain, dummy class */
      MeshDomainDef meshDom = new MeshDomain();
                
      MeshType meshType = new HNMeshType3D();
      MeshSource mesher = new HNMesher3D(0.0, 0.0, 0.0 , 1.2 , 1.2 , 1.2 , 3 , 3 , 3 , meshType , refGeom , meshDom );
      Mesh mesh = mesher.getMesh();
      curve.setMesh(mesh);

      cout << "Nr Points  "<<mesh.numCells(0) << endl;
      cout << "My Rank is :" << myrank << endl;
   
    //WatchFlag watchMe("watch me"); 
    //watchMe.setParam("symbolic preprocessing", 0);
    //watchMe.setParam("discrete function evaluation", 6);
    //watchMe.setParam("integration setup", 6);
    //watchMe.setParam("integral transformation", 6);
    //watchMe.setParam("integration", 6);
    //watchMe.setParam("integration setup", 6); 
    //watchMe.setParam("assembler setup", 6);
    //watchMe.setParam("assembly loop", 6); 
    //watchMe.setParam("matrix config",6);
    //watchMe.setParam("fill", 6);
    //watchMe.setParam("evaluation", 6);
    //watchMe.setParam("dof map setup",6);
    //watchMe.setParam("dof map access", 6); 
      
      // Define the domains 
      CellFilter interior = new MaximalCellFilter(); 
      //CellFilter boundary = new BoundaryCellFilter();
      CellFilter boundary = new DimensionalCellFilter(2);

      // Find left, right, top, bottom boundary edges
      CellFilter left = boundary.subset(new LeftPointTest());
      CellFilter right = boundary.subset(new RightPointTest());
      CellFilter top = boundary.subset(new TopPointTest());
      CellFilter bottom = boundary.subset(new BottomPointTest());
      CellFilter back = boundary.subset(new FrontPointTest());
      CellFilter front = boundary.subset(new BackPointTest());
            
      /* Object needed for curve integral calculation */
      ParametrizedCurve curveIntegral = new ParamCurveIntegral(curve);
      
      CellPredicate curveIN = new CellCurvePredicate( curve , Inside_Curve );
      CellPredicate curveOUT = new CellCurvePredicate( curve , Outside_Curve );
      CellPredicate curveON = new CellCurvePredicate( curve , On_Curve );
      
      //define the different domains 
      CellFilter InCircle = interior.subset(curveIN);
      CellFilter OutsideCircle = interior.subset(curveOUT);
      CellFilter OnCircle = interior.subset(curveON);

      // Define the unknown function and the testfunction 
      Expr vs = new TestFunction(new Lagrange(1),"v");
      Expr us = new UnknownFunction(new Lagrange(1),"u");

      // Define the derivatives
      Expr dx = new Derivative(0); 
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx,dy,dz);
      Expr x1 = new CoordExpr(0);
      Expr x2 = new CoordExpr(1);
      Expr x3 = new CoordExpr(2);

      double C = 50.0;
      
      //double cDirich = 1.0;
      Expr cDirich = new UserDefOp(List(x1,x2,x3), rcp(new CurveExpr(curve,1)) );
      
      Expr h = new CellDiameterExpr();
      Expr nx = new CurveNormExpr(0);
      Expr ny = new CurveNormExpr(1);
      Expr nz = new CurveNormExpr(2);
      Expr alpha = C / h;
      
      // Define the quadrature
      QuadratureFamily quad2 = new GaussianQuadrature(4);
      QuadratureFamily quad_curve_2 = new GaussianQuadrature(6);
      QuadratureFamily quad_curve = new SurfQuadrature(quad_curve_2);
      //QuadratureFamily quad_hi = new GaussianQuadrature(4);
      QuadratureFamily quad_hi = new GaussLobattoQuadrature(2);
      QuadratureFamily tquad = new TrapesoidQuadrature(4);

      // Define the integral 
      Expr eqn = Integral(OutsideCircle, (grad*us)*(grad*vs) - 3*vs, quad2 )
               //+ Integral(InCircle , 1e+3*(1.0-us)*(vs) , quad2 )
	       + Integral(InCircle , (grad*us)*(grad*vs) , quad2 )
               //+ Integral(OnCircle , (grad*us)*(grad*vs) - 3.0*vs , tquad , curve );
	       + Integral(OnCircle , (grad*us)*(grad*vs) - 3.0*vs , quad_hi , curve )
             + Integral(OnCircle , alpha*us*vs - cDirich*alpha*vs, quad_curve , curveIntegral )
- Integral(OnCircle , (nx*(dx*us)+ny*(dy*us)+nz*(dz*us))*vs + us*(nx*(dx*vs)+ny*(dy*vs)+nz*(dz*vs)) , quad_curve , curveIntegral )
          + Integral(OnCircle , cDirich*(nx*(dx*vs)+ny*(dy*vs)+nz*(dz*vs)), quad_curve , curveIntegral );

      // Define the boundary conditions
      Expr bc;

      // We can now set up the linear problem!
      LinearProblem prob(mesh, eqn, bc, vs, us, vecType);

      // Read the parameters for the linear solver from an XML file
#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/bicgstab.xml"));
#else
      ParameterXMLFileReader reader("bicgstab.xml");
#endif
      ParameterList solverParams = reader.getParameters();
   
      // Now we define the problem
      LinearSolver<double> linSolver 
        = LinearSolverBuilder::createSolver(solverParams);
 
      // solve the problem
      Expr soln = prob.solve(linSolver);
      // Write the field in VTK format 

      FieldWriter w = new VTKWriter("Poisson_3D_Surf");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(soln[0])); 
      w.write(); 
      
      Expr err = cDirich - soln[0];
      Expr errExpr = Integral( OnCircle , err*err , quad_curve , curveIntegral );
      FunctionalEvaluator errInt( mesh , errExpr);
      double errorSq = errInt.evaluate();
      Out::os() << "Curve Integration error = " << errorSq << "  , sqrt(error)=" << sqrt(errorSq) <<std::endl;
      Sundance::passFailTest( sqrt(errorSq) , 5e-3);
      
      Expr testExpr = Integral( OnCircle , 1 , quad_curve , curveIntegral );
      FunctionalEvaluator testInt( mesh , testExpr);
      double area = testInt.evaluate();
      Out::os() << " Area of the intersection = " << area <<std::endl;
      Sundance::passFailTest( area - 1.985555 , 1e-5);
      
      Expr testExpr1 = Integral( OutsideCircle , 1 , quad2 )
                       + Integral( OnCircle , 1 , quad_hi , curve );
      FunctionalEvaluator testInt1( mesh , testExpr1);
      double vol = testInt1.evaluate();
      Out::os() << " volume of the intersection = " << vol <<std::endl;
      Sundance::passFailTest( vol - 0.82134 , 1e-5);
      
  }
  catch(std::exception& e)
  {
      Sundance::handleException(e);
  } 
  Sundance::finalize();
}
