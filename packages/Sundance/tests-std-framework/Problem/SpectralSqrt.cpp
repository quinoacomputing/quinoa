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


#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_Playa_Group.hpp"


/*
 * Solve the equation u^2 - alpha = 0, where alpha is a stochastic
 * variable.
 */

bool SpectralSqrt()
{
  int np = MPIComm::world().getNProc();

  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Create a mesh. It will be of type BasisSimplicialMesh, and will
   * be built using a PartitionedLineMesher. */
  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 1*np, meshType);
  Mesh mesh = mesher.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();

  /* Create the Spectral Basis */
  int ndim = 1;
  int order = 2;
  SpectralBasis sbasis = new HermiteSpectralBasis(ndim, order); 
      
  /* Create unknown and test functions, discretized using first-order
   * Lagrange interpolants */
  Expr u = new UnknownFunction(new Lagrange(1), sbasis, "u");
  Expr v = new TestFunction(new Lagrange(1), sbasis, "v");

  /* Create the stochastic input function. */
  Expr a0 = new Sundance::Parameter(1.0);
  Expr a1 = new Sundance::Parameter(0.1);
  Expr a2 = new Sundance::Parameter(0.01);
  Expr alpha = new SpectralExpr(sbasis, tuple(a0, a1, a2));

  /* Create a discrete space, and discretize the function 1.0 on it */
  cout << "forming discrete space" << std::endl;
  DiscreteSpace discSpace(mesh, new Lagrange(1), sbasis, vecType);
  cout << "forming discrete func" << std::endl;
  Expr u0 = new DiscreteFunction(discSpace, 0.5, "u0");

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad = new GaussianQuadrature(2);

  /* Now we set up the weak form of our equation. */
  Expr eqn = Integral(interior, v*(u*u-alpha), quad);

  cout << "equation = " << eqn << std::endl;

  /* There are no boundary conditions for this problem, so the
   * BC expression is empty */
  Expr bc;

  /* We can now set up the nonlinear problem! */

  NonlinearProblem prob(mesh, eqn, bc, v, u, u0, vecType);



#ifdef HAVE_CONFIG_H
  ParameterXMLFileReader reader(searchForFile("SolverParameters/nox.xml"));
#else
  ParameterXMLFileReader reader("nox.xml");
#endif
  ParameterList noxParams = reader.getParameters();

  std::cerr << "solver params = " << noxParams << std::endl;

  NOXSolver solver(noxParams);

  prob.solve(solver);

  /* Inspect solution values. The solution is constant in space,
   * so we can simply take the first NTerms entries in the vector */
  Vector<double> vec = DiscreteFunction::discFunc(u0)->getVector();

  for (int i=0; i<vec.space().numLocalElements(); i++)
  {
    cout << "u[" << i << "] = " << vec[i] << std::endl;
  }

  double tol = 1.0e-12;
  double errorSq = 0.0;
  return SundanceGlobal::checkTest(errorSq, tol);
}
