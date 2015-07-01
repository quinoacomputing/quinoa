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
 * Projects a function onto a high-order basis
 */


bool HighOrderProjection2D()
{
  int np = MPIComm::world().getNProc();

  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Create a mesh. It will be of type BasisSimplicialMesh, and will
   * be built using a PartitionedRectangleMesher. */
  int nx = 4;
  int ny = 2;
  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, np,
    0.0, 1.0, ny, 1,
    meshType);
  Mesh mesh = mesher.getMesh();


  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
  CellFilter edges = new DimensionalCellFilter(1);

  CellFilter left = edges.subset(new CoordinateValueCellPredicate(0,0.0));
  CellFilter right = edges.subset(new CoordinateValueCellPredicate(0,1.0));
  CellFilter top = edges.subset(new CoordinateValueCellPredicate(1,1.0));
  CellFilter bottom = edges.subset(new CoordinateValueCellPredicate(1,0.0));

      
  /* Create unknown and test functions, discretized using first-order
   * Lagrange interpolants */
  int order = 2;
  double p = (double) order;
  Expr u = new UnknownFunction(new Lagrange(order), "u");
  Expr v = new TestFunction(new Lagrange(order), "v");

  /* Create coordinate functions */
  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad6 = new GaussianQuadrature(6);

  /* Define the weak form */
  Expr alpha = sqrt(1.0);
  Expr z = x + alpha*y;
  Expr exactSoln = pow(z, p);
  Expr eqn = Integral(interior, (v*(u - exactSoln)), quad6);

  /* Define the Dirichlet BC */
  Expr bc;

  /* We can now set up the linear problem! */
  LinearProblem prob(mesh, eqn, bc, v, u, vecType);

  LinearSolver<double> solver 
    = LinearSolverBuilder::createSolver("aztec.xml");

  Expr soln = prob.solve(solver);

  /* compute the error */
  Expr err = exactSoln - soln;
  Expr errExpr = Integral(interior, 
    err*err,
    quad6);


  FunctionalEvaluator errInt(mesh, errExpr);

  double errorSq = errInt.evaluate();
  cout << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

  return SundanceGlobal::checkTest(sqrt(errorSq), 1.0e-11);
}

