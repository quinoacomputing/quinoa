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
#include "SundanceUnknownParameter.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionData.hpp"


/** 
 * Solves the steady Burgers equation in 1D with forcing
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
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 100*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter rightPoint = points.subset(new RightPointTest());
      CellFilter leftPoint = points.subset(new LeftPointTest());
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* Create a discrete space, and discretize the function 1.0 on it */
      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      L2Projector projector(discSpace, 1.0);
      Expr u0 = projector.project();

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(4);

      /* Parameters */
      Expr p = new UnknownParameter("p");
      Expr p0 = new Sundance::Parameter(2.0);

      /* Forcing term */
      Expr f = p * (p*x*(2.0*x*x - 3.0*x + 1.0) + 2.0);

      /* Define the weak form */
      Expr eqn = Integral(interior, (dx*u)*(dx*v) + v*u*(dx*u) - v*f, quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint+rightPoint, v*u, quad);

      /* Create a Playa NonlinearOperator object */
      NonlinearProblem prob(mesh, eqn, bc, v, u, u0, 
        p, p0, vecType);

      ParameterXMLFileReader reader("nox.xml");
      ParameterList noxParams = reader.getParameters();

      NOXSolver solver(noxParams);
      LinearSolver<double> linSolver = solver.linSolver();

      // Solve the nonlinear system
      SolverState<double> status = prob.solve(solver);
      TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
        runtime_error, "solve failed");


      /* compute senstivities */
      Expr sens = prob.computeSensitivities(linSolver);

      /* Write the field in ASCII format */
      FieldWriter w = new MatlabWriter("Burgers1DSoln");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(u0));
      w.addField("sens_a", new ExprFieldWrapper(sens[0]));
      w.write();


      /* check solution */
      Expr errExpr = Integral(interior, 
        pow(u0-p0*x*(1.0-x), 2),
        new GaussianQuadrature(8));
      Expr errExprA = Integral(interior, 
        pow(sens[0]-x*(1.0-x), 2),
        new GaussianQuadrature(8));

      double errorSq0 = evaluateIntegral(mesh, errExpr);
      std::cerr << "soln error norm = " << sqrt(errorSq0) << std::endl << std::endl;

      double errorSqA = evaluateIntegral(mesh, errExprA);
      std::cerr << "sens A error norm = " << sqrt(errorSqA) << std::endl << std::endl;

      double error = sqrt(errorSq0 + errorSqA);
      
      double tol = 1.0e-4;
      Sundance::passFailTest(error, tol);

    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
