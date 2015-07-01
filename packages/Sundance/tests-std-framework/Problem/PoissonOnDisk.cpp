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
#include "SundanceCellVectorExpr.hpp"
#include "SundanceEvaluator.hpp"

using Sundance::List;
/** 
 * Solves the Poisson equation in 2D on the unit disk
 */


bool PoissonOnDisk()
{
#ifdef HAVE_SUNDANCE_EXODUS

  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Get a mesh */
  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource meshReader = new ExodusNetCDFMeshReader("disk.ncdf", meshType);
  Mesh mesh = meshReader.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
  CellFilter bdry = new BoundaryCellFilter();


      
  /* Create unknown and test functions, discretized using first-order
   * Lagrange interpolants */
  Expr u = new UnknownFunction(new Lagrange(1), "u");
  Expr v = new TestFunction(new Lagrange(1), "v");

  /* Create differential operator and coordinate functions */
  Expr dx = new Derivative(0);
  Expr dy = new Derivative(1);
  Expr grad = List(dx, dy);
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad2 = new GaussianQuadrature(2);
  QuadratureFamily quad4 = new GaussianQuadrature(4);

  /* Define the weak form */
  Expr eqn = Integral(interior, (grad*v)*(grad*u)  + v, quad2);
  /* Define the Dirichlet BC */
  Expr bc = EssentialBC(bdry, v*u, quad4);


  /* We can now set up the linear problem! */
  LinearProblem prob(mesh, eqn, bc, v, u, vecType);

  LinearSolver<double> solver 
    = LinearSolverBuilder::createSolver("amesos.xml");

  Expr soln = prob.solve(solver);

  double R = 1.0;
  Expr exactSoln = 0.25*(x*x + y*y - R*R);

  DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
  Expr du = L2Projector(discSpace, exactSoln-soln).project();

  /* Write the field in VTK format */
  FieldWriter w = new VTKWriter("PoissonOnDisk");
  w.addMesh(mesh);
  w.addField("soln", new ExprFieldWrapper(soln[0]));
  w.addField("error", new ExprFieldWrapper(du));
  w.write();

      
  Expr errExpr = Integral(interior, 
    pow(soln-exactSoln, 2),
    new GaussianQuadrature(4));

  double errorSq = evaluateIntegral(mesh, errExpr);
  std::cerr << "soln error norm = " << sqrt(errorSq) << std::endl << std::endl;


  /* Check error in automatically-computed cell normals */
  /* Create a cell normal expression. Note that this is NOT a constructor
   * call, hence no "new" before the CellNormalExpr() function. The
   * argument "2" is the spatial dimension (mandatory), and 
   * the "n" is the name of the expression (optional). 
   */
  Expr n = CellNormalExpr(2, "n");
  Expr nExact = List(x, y)/sqrt(x*x + y*y);
  Expr nErrExpr = Integral(bdry, pow(n-nExact, 2.0), new GaussianQuadrature(1));
  double nErrorSq = evaluateIntegral(mesh, nErrExpr);
  Out::root() << "normalVector error norm = " 
              << sqrt(nErrorSq) << std::endl << std::endl;

  double tol = 1.0e-4;
  return SundanceGlobal::checkTest(sqrt(errorSq + nErrorSq), tol);
#else
  std::cout << "dummy PoissonOnDisk PASSED. Enable Exodus to run the actual test" << std::endl;
  return true;
#endif
}
