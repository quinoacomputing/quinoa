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

/** \example Scheibe.cpp
* Solves the elastic panel example with 2 elements and symmetric
condition from FE section of Ramm's Baustatik script. */
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(BottomPoint, {return (fabs(x[1]) < 1.0e-10)&&(fabs(x[0]-30.0) < 1.0);})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-60.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-80.0) < 1.0e-10;})
CELL_PREDICATE(TopPointCheckTest, {return (fabs(x[0]-30.0) < 2.0)&&(fabs(x[1]-80.0) < 2.0);} )
int main(int argc, char** argv)
{
try
{
        Sundance::init(&argc, &argv);
        /* Linear algebra with Epetra */
        VectorType<double> vecType = new EpetraVectorType();
	
        /* Create one circle */
        ParametrizedCurve circle = new Circle( 30.0 , 40.0 , 20.0 , 1 , 0.00001 );
	
        /* Meshing with PartitionedRectangleMesher. Gives only simplices
           unfortunately. We work on a nx1 x nx2 grid */
        int nx1 = 15;
        int nx2 = 20;
        MeshType meshType = new BasicSimplicialMeshType();
        MeshSource mesher = new PartitionedRectangleMesher(0.0,60.0,nx1,1,0.0,80.0,nx2,1,meshType);
        Mesh mesh = mesher.getMesh();
	
    WatchFlag watchMe("watch me");
    watchMe.setParam("symbolic preprocessing", 0);
    watchMe.setParam("discrete function evaluation", 0);
    watchMe.setParam("integration setup", 1);
    watchMe.setParam("integral transformation", 1);
    watchMe.setParam("integration", 6);
    watchMe.setParam("fill", 0);
    watchMe.setParam("evaluation", 0);
	
        /* Cell filter that contains the complete domain */
        CellFilter Omega = new MaximalCellFilter();
	
        /* Cell Filter that contains the boundaries */
        CellFilter edges = new DimensionalCellFilter(1);
        CellFilter points = new DimensionalCellFilter(0);

        /* Find free boundary edges */
        CellFilter left = edges.subset(new LeftPointTest());
        CellFilter symmetry = edges.subset(new RightPointTest());
        CellFilter top = edges.subset(new TopPointTest());
        CellFilter bottom = edges.subset(new BottomPointTest());
	  
        /* Displacement based shape functions + test functions
        (Linear Lagrangian interpolants -> CST) */
        BasisFamily La1 = new Lagrange(1);
        Expr ux = new UnknownFunction(La1);
        Expr uy = new UnknownFunction(La1);
        Expr dux = new TestFunction(La1);
        Expr duy = new TestFunction(La1);
	
        Expr u = List(ux , uy );
        Expr du = List( dux , duy );

        /* Create differential operator */
        Expr x1 = new CoordExpr(0);
        Expr x2 = new CoordExpr(1);
        Expr dx1 = new Derivative(0);
        Expr dx2 = new Derivative(1);
        Expr L1 = List(dx1,0.0);
        Expr L2 = List(0.0,dx2);
        Expr L3 = List(dx2,dx1);
        Expr L = List(L1,L2,L3);
	
        CellPredicate curveIN = new CellCurvePredicate( circle , Inside_Curve );
        CellPredicate curveOUT = new CellCurvePredicate( circle , Outside_Curve );
        CellPredicate curveON = new CellCurvePredicate( circle , On_Curve );
      
        CellFilter InCircle = Omega.subset(curveIN);
        CellFilter OutsideCircle = Omega.subset(curveOUT);
        CellFilter OnCircle = Omega.subset(curveON);
      
        /* Plane Stress Constitutive Relations */
        // 200 6000 900 ,
        // nu 0.29
        // dicke = 1
        // F = 100
        Expr nu = new Sundance::Parameter(0.29, "nu");
        Expr E = new Sundance::Parameter(206900, "E");
        Expr Dicke = new Sundance::Parameter(1, "Dicke");
	
        // we just add one expression in order to trigger QuadratureIntegral to happen
        // so that we can test QuadratureIntegral as well for the ACI logic
        // If we add this to the E, then we should get almost the same solution
        Expr coordExpr = new Sundance::CoordExpr(0);

        Expr C = ((coordExpr + E*Dicke)/(1-nu*nu))*List(List(1.0, nu , 0.0),List( nu , 1.0, 0.0),List(0.0, 0.0, 0.5-0.5*nu));
        Expr C_FCM = 0.00001*C;
        double  F = -100;
	
        /* Quadrature rule for numerical integration (not needed here)*/
        QuadratureFamily quad_hi = new GaussianQuadrature(10);
        QuadratureFamily quad = new GaussianQuadrature(2);
	
        /* Weak form of elastostatic problem (Principle of Virtual Displacements) */
        Expr Ldu = L*du;
        Expr LduT= List(Ldu);
        Expr eqn = Integral(OutsideCircle, LduT*C*(L*u), quad ) +
	               Integral(OnCircle, LduT*C*(L*u), quad_hi , circle ) +
		           Integral(InCircle, LduT*C_FCM*(L*u), quad ) +
	               Integral(top,F*du[1], quad );
	
        /* Define the Dirichlet BC */
        Expr bc = EssentialBC(bottom, du[1]*u[1], quad) +
		          EssentialBC(bottom, du[0]*u[0], quad);
	
        /* We can now set up the linear problem! */
        LinearProblem prob(mesh, eqn, bc, du, u, vecType);
	
        /* Read the parameters for the linear solver from an XML file
        This file can be found under examples-tutorials in the sundance
        directory*/
        ParameterXMLFileReader reader("bicgstab.xml");
        ParameterList solverParams = reader.getParameters();
        LinearSolver<double> linSolver = LinearSolverBuilder::createSolver(solverParams);

        /* solve the problem */
        Expr soln = prob.solve(linSolver);
	
        /* Project the stresses onto a discrete space */
        DiscreteSpace discreteSpace(mesh, List( La1 , La1 , La1 ), vecType);
        L2Projector projector(discreteSpace, C*(L*soln));
        Expr stresses = projector.project();
	
        /* Write the field in VTK format */
        FieldWriter w = new VTKWriter("PerforatedPlate");
        w.addMesh(mesh);
        w.addField("u", new ExprFieldWrapper(soln[0]));
        w.addField("v", new ExprFieldWrapper(soln[1]));
        w.addField("sigxx", new ExprFieldWrapper(stresses[0]));
        w.addField("sigyy", new ExprFieldWrapper(stresses[1]));
        w.addField("tauxy", new ExprFieldWrapper(stresses[2]));
        w.write();

        /* ========== Benchmark values testing ============= */
        CellFilter checkPoint = points.subset(new TopPointCheckTest());

        Expr checkVerticalDispErr = Integral(checkPoint, (soln[1]-0.121)*(soln[1]-0.121),
        		new GaussianQuadrature(2) );
        Expr checkHorisontalDispErr = Integral(checkPoint, (soln[0])*(soln[0]) ,
        		new GaussianQuadrature(2) );

        FunctionalEvaluator verticalDisp(mesh, checkVerticalDispErr);
        FunctionalEvaluator horizontalDisp(mesh, checkHorisontalDispErr);

        double VerticalDispErr = verticalDisp.evaluate();
        double HorisontalDispErr = horizontalDisp.evaluate();

        std::cout << "VerticalDispErr = " << VerticalDispErr << std::endl;
        std::cout << "HorisontalDispErr = " << HorisontalDispErr << std::endl;

        double tol = 5.0e-2;
        Sundance::passFailTest(sqrt(VerticalDispErr)+sqrt(HorisontalDispErr), tol);
}
catch(std::exception& e)
{
      Sundance::handleException(e);
}
Sundance::finalize();
}
