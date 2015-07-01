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
#include "SundancePeriodicLineMesher.hpp"
#include "SundancePeriodicMeshType1D.hpp"

/** 
 * Solves the equation
 * 
 * \f[ u'' + 2 u' + u^2 = \sin^2(2x+1) - 4 \sin(2x+1) + 4 \cos(2x+1) \f]
 *
 * with periodic BC on the interval \f$ (0, 2\pi) \f$.
 */

bool NonlinearPeriodic1D()
{
  int np = MPIComm::world().getNProc();
  TEUCHOS_TEST_FOR_EXCEPT(np != 1);

  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Create a periodic mesh */
  int nx = 1000;
  const double pi = 4.0*atan(1.0);
  MeshType meshType = new PeriodicMeshType1D();
  MeshSource mesher = new PeriodicLineMesher(0.0, 2.0*pi, nx, meshType);
  Mesh mesh = mesher.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
      
  /* Create unknown and test functions, discretized using first-order
   * Lagrange interpolants */
  Expr u = new UnknownFunction(new Lagrange(1), "u");
  Expr v = new TestFunction(new Lagrange(1), "v");

  /* Create differential operator and coordinate function */
  Expr dx = new Derivative(0);
  Expr x = new CoordExpr(0);

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad = new GaussianQuadrature(4);

      
  /* Define the weak form */
  Expr s = sin(2.0*x+1.0);
  Expr c = cos(2.0*x+1.0);
  Expr eqn = Integral(interior, 
    (dx*v)*(dx*u) - 2.0*v*(dx*u) - v*u*u + v*(s*s - 4.0*s + 4.0*c),
    quad);
  Expr bc ; // no explicit BC needed

  /* We can now set up the linear problem! */

  DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
  Expr u0 = new DiscreteFunction(discSpace, 1.0);

  NonlinearProblem prob(mesh, eqn, bc, v, u, u0, vecType);


  ParameterXMLFileReader reader("nox.xml");
  ParameterList solverParams = reader.getParameters();

  NOXSolver solver(solverParams);
  prob.solve(solver);

  Expr uExact = s;

  Expr uErr = uExact - u0;
      
  Expr uErrExpr = Integral(interior, 
    uErr*uErr,
    new GaussianQuadrature(6));
      
  FunctionalEvaluator uErrInt(mesh, uErrExpr);

  double uErrorSq = uErrInt.evaluate();
  std::cerr << "u error norm = " << sqrt(uErrorSq) << std::endl << std::endl;

  double tol = 1.0e-3;
  return SundanceGlobal::checkTest(sqrt(uErrorSq), tol);
}
