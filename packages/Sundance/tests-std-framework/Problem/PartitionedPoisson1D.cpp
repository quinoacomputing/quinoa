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


/** 
 * Solves the Poisson equation in 1D
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})

  int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);

    int np = MPIComm::world().getNProc();

    if (np > 1)
    {
      cout << "parallel partioned poisson test INACTIVE" << std::endl;
      return 0;
    }
    else
    {

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      //      PartitionedLineMesher::classVerbosity() = VerbExtreme;

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 10*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2);


      /* Define the weak form */
      Expr eqn = Integral(interior, -(dx*v)*(dx*u) - 2.0*v, quad);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, v*u, quad);



      /* We can now set up the linear problem! */
      ParameterList vp = *LinearProblem::defaultVerbParams();
      LinearProblem prob(mesh, eqn, bc, v, u, vecType, vp, true); 

      LinearOperator<double> A = prob.getOperator();
      Vector<double> b = prob.getRHS()[0];
      b.scale(-1.0);

      cout << "b=" << b << std::endl;

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/aztec-ifpack.xml"));;
#else
      ParameterXMLFileReader reader("aztec-ifpack.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      cout << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Vector<double> bi = b.getBlock(0);
      Vector<double> bb = b.getBlock(1);
      LinearOperator<double> Aii = A.getBlock(0,0);
      LinearOperator<double> Aib = A.getBlock(0,1);
      LinearOperator<double> Abb = A.getBlock(1,1);
      Vector<double> xi = bi.copy();
      Vector<double> xb = bb.copy();

      solver.solve(Abb, bb, xb);
      Vector<double> c = bi - Aib * xb;
      solver.solve(Aii, c, xi);

      cout << "x_i=" << std::endl << xi << std::endl;
      cout << "x_b=" << std::endl << xb << std::endl;

      Vector<double> solnVec = prob.convertToMonolithicVector(tuple(xi), tuple(xb));
      Expr soln = prob.formSolutionExpr(tuple(solnVec));

      Expr exactSoln = x*(x-2.0);

      Expr errExpr = Integral(interior, 
        pow(soln-exactSoln, 2),
        new GaussianQuadrature(4));

      Expr derivErr = dx*(exactSoln-soln);
      Expr derivErrExpr = Integral(interior, 
        pow(dx*(soln-exactSoln), 2),
        new GaussianQuadrature(2));

      Expr fluxErrExpr = Integral(leftPoint, 
        pow(dx*(soln-exactSoln), 2),
        new GaussianQuadrature(2));

      double errorSq = evaluateIntegral(mesh, errExpr);
      cout << "error norm = " << sqrt(errorSq) << std::endl << std::endl;



      double derivErrorSq = evaluateIntegral(mesh, derivErrExpr);
      cout << "deriv error norm = " << sqrt(derivErrorSq) << std::endl << std::endl;

      double fluxErrorSq = evaluateIntegral(mesh, fluxErrExpr);
      cout << "flux error norm = " << sqrt(fluxErrorSq) << std::endl << std::endl;

      Expr exactFluxExpr = Integral(leftPoint, 
        dx*exactSoln,
        new GaussianQuadrature(2));

      Expr numFluxExpr = Integral(leftPoint, 
        dx*soln,
        new GaussianQuadrature(2));

      cout << "exact flux = " << evaluateIntegral(mesh, exactFluxExpr) << std::endl;
      cout << "computed flux = " << evaluateIntegral(mesh, numFluxExpr) << std::endl;


      double tol = 1.0e-12;
      Sundance::passFailTest(sqrt(errorSq + derivErrorSq + fluxErrorSq), tol);
    }
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}
