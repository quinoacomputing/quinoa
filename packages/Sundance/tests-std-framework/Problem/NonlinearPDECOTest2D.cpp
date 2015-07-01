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

#include "SundanceZeroExpr.hpp"

CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;}) 
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;}) 
CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;}) 
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;}) 


int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();
      int npx = -1;
      int npy = -1;
      PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
      TEUCHOS_TEST_FOR_EXCEPT(npx < 1);
      TEUCHOS_TEST_FOR_EXCEPT(npy < 1);
      TEUCHOS_TEST_FOR_EXCEPT(npx * npy != np);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 32;
      int ny = 32;
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, npx, 
        0.0,  1.0, ny, npy, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter(); 
      CellFilter edges = new DimensionalCellFilter(1);
      
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());
      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());

      /* Create unknown function */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr lambda = new UnknownFunction(new Lagrange(1), "lambda");
      Expr alpha = new UnknownFunction(new Lagrange(1), "alpha");
      
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      
      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      WatchFlag watch("watch me");
      watch.setParam("integration setup", 6);
      watch.deactivate();

      const double pi = 4.0*atan(1.0);

      double R = 0.001;

      /* uStar, spliced in from Mathematica */
      Expr uStar = 2*pow(pi,2)*R*pow(cos(pi*y),2)*
        pow(sin(pi*x),2) + 
        sin(pi*y)*(sin(pi*x) + 
          2*pow(pi,2)*R*
          pow(cos(pi*x),2)*sin(pi*y)
          - 2*pow(pi,2)*R*
          pow(sin(pi*x),2)*sin(pi*y)
          + 2*R*pow(sin(pi*x),3)*
          pow(sin(pi*y),2));
      
      Expr mismatch = u-uStar;
      Expr fit = Integral(interior, 0.5*mismatch*mismatch, quad, watch);
      Expr reg = Integral(interior, 0.5*R*alpha*alpha, quad, watch);

      Expr g = 2.0*pi*pi*u + u*u;

      Expr constraint = Integral(interior, (grad*u)*(grad*lambda) - lambda*g + lambda*alpha, quad, watch);

      Expr lagrangian = fit + reg + constraint;

      Expr bc = EssentialBC(top+bottom+left+right, 
        lambda*u, quad, watch);

      BasisFamily L1=new Lagrange(1);
      DiscreteSpace discSpace(mesh, List(L1, L1, L1), vecType);
      Expr W0 = new DiscreteFunction(discSpace, 0.0);
      Expr u0 = W0[0];
      Expr lambda0 = W0[1];
      Expr alpha0 = W0[2];
      Expr W = List(u, lambda, alpha);

      Expr dummy;

      Functional L(mesh, lagrangian, bc, vecType);
      
      NonlinearProblem prob 
        = L.nonlinearVariationalProb(W, W0, W, W0, dummy, dummy);

      ParameterXMLFileReader reader("nox-amesos.xml");

      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams);
      
      SolverState<double> status = prob.solve(solver);
      TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
        runtime_error, "solve failed");


      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("pdeco2D");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(W0[0]));
      w.addField("lambda", new ExprFieldWrapper(W0[1]));
      w.addField("alpha", new ExprFieldWrapper(W0[2]));
      w.write();

      Expr uExact = sin(pi*x)*sin(pi*y);
      Expr alphaExact = uExact*uExact;
      Expr lambdaExact = -R*alphaExact;

      Expr uErrExpr = Integral(interior, 
        pow(u0-uExact, 2),
        new GaussianQuadrature(8));

      Expr alphaErrExpr = Integral(interior, 
        pow(alpha0-alphaExact, 2),
        new GaussianQuadrature(8));

      Expr lambdaErrExpr = Integral(interior, 
        pow(lambda0-lambdaExact, 2),
        new GaussianQuadrature(8));

      
      double uErrorSq = evaluateIntegral(mesh, uErrExpr);
      std::cerr << "u error norm = " << sqrt(uErrorSq) << std::endl << std::endl;

      double alphaErrorSq = evaluateIntegral(mesh, alphaErrExpr);
      std::cerr << "alpha error norm = " << sqrt(alphaErrorSq) << std::endl << std::endl;

      double lambdaErrorSq = evaluateIntegral(mesh, lambdaErrExpr);
      std::cerr << "lambda error norm = " << sqrt(lambdaErrorSq) << std::endl << std::endl;

      double err = sqrt(uErrorSq + lambdaErrorSq + alphaErrorSq);

      double tol = 0.05;
      Sundance::passFailTest(err, tol);

    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize();
  return Sundance::testStatus(); 
}
