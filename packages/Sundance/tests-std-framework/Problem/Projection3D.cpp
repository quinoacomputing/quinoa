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
#include "SundanceEvaluator.hpp"

using Sundance::List;
/** 
 * Projects a function onto a high-order basis
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})


#if defined(HAVE_SUNDANCE_EXODUS) && defined(Trilinos_DATA_DIR)


int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher 
        = new ExodusMeshReader("cube-0.1", meshType);
      Mesh mesh = mesher.getMesh();


      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      int order = 2;
      double p = (double) order;
      Expr u = new UnknownFunction(new Lagrange(order), "u");
      Expr v = new TestFunction(new Lagrange(order), "v");

      /* Create coordinate functions */
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad6 = new GaussianQuadrature(6);

      /* Define the weak form */
      Expr alpha = sqrt(2.0);
      Expr beta = sqrt(3.0);
      Expr w = x + alpha*y + beta*z;
      Expr exactSoln = pow(w, p);
      Expr eqn = Integral(interior, (v*(u - exactSoln)), quad6);

      /* Define the Dirichlet BC */
      Expr bc;

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);


      ParameterXMLFileReader reader("amesos.xml");
      ParameterList solverParams = reader.getParameters();

      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver);

      /* compute the error */
      Expr err = exactSoln - soln;
      Expr errExpr = Integral(interior, 
                              err*err,
                              quad6);


      FunctionalEvaluator errInt(mesh, errExpr);

      double errorSq = errInt.evaluate();
      cout << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

      Sundance::passFailTest(sqrt(errorSq), 1.0e-11);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}

#else

int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy Projection3D PASSED. Enable exodus to run the actual test" <<
 std::endl;
  Sundance::finalize();
  return 0;
}

#endif
