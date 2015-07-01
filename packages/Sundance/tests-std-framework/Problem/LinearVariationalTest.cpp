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

using Sundance::List;
/** 
 *
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 40;
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);

      CellFilter right = points.subset(new RightPointTest());
      CellFilter left = points.subset(new LeftPointTest());


      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr alpha = new UnknownFunction(new Lagrange(2), "alpha");
      Expr lambda = new UnknownFunction(new Lagrange(2), "lambda");

      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);

      /* set initial values for u0, alpha0, and lambda0 */
      Expr u0 = new DiscreteFunction(discSpace, 0.0, "u0");

      std::cerr << "doing projection " << std::endl;
      L2Projector projector(discSpace, x*(1.0-x));
      Expr alpha0 = projector.project();

      Expr lambda0 = new DiscreteFunction(discSpace, 0.0, "lambda0");

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily q2 = new GaussianQuadrature(2);
      QuadratureFamily q4 = new GaussianQuadrature(4);

      /* the target we're trying to match */
      const double pi = 4.0*atan(1.0);
      Expr target = sin(pi*x);

      /* regularization weight */
      double R = 0.1;

      /* Define the Lagrangian */
      Expr lag = Integral(interior, 0.5*pow(u - target, 2.0) + R*alpha*alpha
                          + (dx*lambda)*(dx*u) + lambda*alpha, q2);

      /* Define the Dirichlet BC */
      Expr lagBC = EssentialBC(left+right, lambda*u, q2); 

      Functional L(mesh, lag, lagBC, vecType);


      /* derive and solve the state equation */

      LinearProblem stateProb = L.linearVariationalProb(lambda, lambda0,
                                                        u,
                                                        alpha, alpha0);

      


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/bicgstab.xml"));
#else
      ParameterXMLFileReader reader("bicgstab.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      std::cerr << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      std::cerr << "solving state equation" << std::endl;

      Expr uSoln = stateProb.solve(solver);

      DiscreteFunction* df = DiscreteFunction::discFunc(uSoln);
      Vector<double> uVec = df->getVector();
      DiscreteFunction::discFunc(u0)->setVector(uVec);

      Expr uExact = -1.0/12.0 * x * (1.0 - 2.0*x*x + x*x*x);
      L2Projector p2(discSpace, uExact);
      Expr discExact = p2.project();
      

      


      Expr uErrExpr = Integral(interior, 
                              pow(uSoln-uExact, 2),
                              new GaussianQuadrature(8));

      double uErrorSq = evaluateIntegral(mesh, uErrExpr);
      std::cerr << "state error norm = " << sqrt(uErrorSq) << std::endl << std::endl;

      /* derive and solve the adjoint equation */
      
      LinearProblem adjointProb = L.linearVariationalProb(u, u0,
                                                          lambda,
                                                          alpha, alpha0);
      std::cerr << "solving adjoint equation" << std::endl;

      Expr lambdaSoln = adjointProb.solve(solver);

      
      Expr lambdaExact 
        = -1.0/12.0 * (pow(x,3)/6.0 - pow(x,5)/10.0 + pow(x,6)/30.0)
        + sin(pi*x)/pi/pi + x/120.0;
      L2Projector p3(discSpace, lambdaExact);
      Expr discLambdaExact = p3.project();



      Expr lambdaErrExpr = Integral(interior, 
                              pow(lambdaSoln-lambdaExact, 2),
                              new GaussianQuadrature(12));

      double lambdaErrorSq = evaluateIntegral(mesh, lambdaErrExpr);
      std::cerr << "adjoint error norm = " << sqrt(lambdaErrorSq) << std::endl << std::endl;

      /* Write the field in Matlab format */
      FieldWriter w = new MatlabWriter("LinearVarTest.dat");
      w.addMesh(mesh);
      w.addField("numerical U", new ExprFieldWrapper(uSoln[0]));
      w.addField("exact U", new ExprFieldWrapper(discExact[0]));
      w.addField("numerical U", new ExprFieldWrapper(lambdaSoln[0]));
      w.addField("exact U", new ExprFieldWrapper(discLambdaExact[0]));
      w.write();
      
      double tol = 1.0e-6;
      
      Sundance::passFailTest(sqrt(lambdaErrorSq + uErrorSq), tol);
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
