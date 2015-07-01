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
   - Laplace u = 0
   stationary
   With Dirichlet boundary condition:
             u = const , at the boundary
 */
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-4;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-4;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-4;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-4;})

REFINE_MESH_ESTIMATE(MeshRefEst , { return 0; } , {return 1;} )
MESH_DOMAIN( MeshDomain , {return true;})

int main(int argc, char** argv)
{
  try
  {
      Sundance::init(&argc, &argv);
      int myrank = MPIComm::world().getRank();

      // We will do our linear algebra using Epetra 
      VectorType<double> vecType = new EpetraVectorType();

    RefinementClass refCl = new MeshRefEst();
    MeshDomainDef meshDom = new MeshDomain();
      
    ParametrizedCurve curve = new Circle( 0.5 , 0.5 , 0.23 , 1 , 1 );
    //ParametrizedCurve curve = new Box2D( 0.15 , 0.15 , 0.10 , 0.10 , 1 , 0.00001 );
      
    /* Object needed for curve integral calculation */
    ParametrizedCurve curveIntegral = new ParamCurveIntegral(curve);
    
    /* the cell predicate for different integrations */
    CellPredicate curveIN = new CellCurvePredicate( curve , Inside_Curve );
    CellPredicate curveOUT = new CellCurvePredicate( curve , Outside_Curve );
    CellPredicate curveON = new CellCurvePredicate( curve , On_Curve );
    
    MeshType meshType = new HNMeshType2D();
    MeshSource mesher = new HNMesher2D(0.0, 0.0, 1.0 , 1.0 , 51 , 51 , meshType , refCl , meshDom );
    Mesh mesh = mesher.getMesh();
   
      // set the Verbosity high for the grid
      //Mesh::classVerbosity() = VerbExtreme; 
    WatchFlag watchMe("watch me");
    watchMe.setParam("symbolic preprocessing", 0);
    watchMe.setParam("discrete function evaluation", 6);
    watchMe.setParam("integration setup", 1);
    watchMe.setParam("integral transformation", 1);
    watchMe.setParam("integration", 6);
    watchMe.setParam("fill", 6);
    watchMe.setParam("evaluation", 6);


      cout << "Nr Points  "<<mesh.numCells(0) << endl;
      cout << "My Rank is :" << myrank << endl;
      
      // Define the domains 
      CellFilter interior = new MaximalCellFilter();
      CellFilter boundary = new BoundaryCellFilter();

      // Find left, right, top, bottom boundary edges
      CellFilter left = boundary.subset(new LeftPointTest());
      CellFilter right = boundary.subset(new RightPointTest());
      CellFilter top = boundary.subset(new TopPointTest());
      CellFilter bottom = boundary.subset(new BottomPointTest());
      
      //define the different domains 
      CellFilter InCircle = interior.subset(curveIN);
      CellFilter OutsideCircle = interior.subset(curveOUT);
      CellFilter OnCircle = interior.subset(curveON);

      
      BasisFamily La = new Lagrange(1);

      // Define the unknown function and the testfunction 
      Expr vs = new TestFunction( La , "v" );
      Expr us = new UnknownFunction( La , "u" );

      //DiscreteSpace MySpace(mesh, La , vecType);
      //Expr us = new DiscreteFunction(MySpace, 0.0, "u");
      
      // Define the derivatives
      Expr dx = new Derivative(0); 
      Expr dy = new Derivative(1);
      Expr grad = List(dx,dy);

      // Define the quadrature
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad_hi = new GaussianQuadrature(10);
      
      double C = 6.0;
      double cDirich = 1.0;
      Expr h = new CellDiameterExpr();
      Expr nx = new CurveNormExpr(0);
      Expr ny = new CurveNormExpr(1);
      Expr alpha = C / h; 
      
      // Define the integral 
      Expr eqn = Integral(OutsideCircle , (grad*us)*(grad*vs) - 3.0*vs , quad2 )
       + Integral(InCircle , 1.0*(grad*us)*(grad*vs) - 1.0*3.0*vs , quad2 )
       + Integral(OnCircle , (grad*us)*(grad*vs) - 3.0*vs , quad_hi , curve )
       // here comes the Nitsche formula for Dirichlet Boundary Condition
       + Integral(OnCircle , alpha*us*vs - cDirich*alpha*vs, quad2 , curveIntegral )
       - Integral(OnCircle , (nx*(dx*us)+ny*(dy*us))*vs + us*(nx*(dx*vs)+ny*(dy*vs)) , quad2 , curveIntegral )
       + Integral(OnCircle , cDirich*(nx*(dx*vs)+ny*(dy*vs)), quad2 , curveIntegral );

      // Define the boundary conditions
      Expr bc =  EssentialBC( right , vs*(us-cDirich) , quad2 );

      // We can now set up the linear problem!
      LinearProblem prob(mesh, eqn , bc, vs, us, vecType);

      // Read the parameters for the linear solver from an XML file

      ParameterXMLFileReader reader("bicgstab.xml");
      //ParameterXMLFileReader reader("superlu_dist.xml");

      ParameterList solverParams = reader.getParameters();
      // Now we define the problem
      LinearSolver<double> linSolver 
        = LinearSolverBuilder::createSolver(solverParams);
 
      // solve the problem
      Expr soln = prob.solve(linSolver);
      // Write the field in VTK format 
      FieldWriter w = new VTKWriter("Poisson_HN_Nitsche");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(soln[0]));
      w.write(); 
      
      // integrate on the curve 
      Evaluator::classVerbosity() = 6;
      Expr err = cDirich - soln;
      Expr errExpr = Integral( OnCircle , err*err , quad2 , curveIntegral , watchMe);
      //Expr errExpr = Integral( left , err*err , quad2 , watchMe );
      FunctionalEvaluator errInt( mesh , errExpr);
      double errorSq = errInt.evaluate();
      Out::os() << "Curve Integration error = " << errorSq << "  , sqrt(error)=" << sqrt(errorSq) <<std::endl;
      
      double tol = 2.5e-1;
      Sundance::passFailTest( sqrt(errorSq) , tol);
  }
  catch(std::exception& e)
  {
      Sundance::handleException(e);
  }
  Sundance::finalize();
}
