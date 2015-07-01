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

#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_Playa_Group.hpp"


#include "PlayaNOXSolver.hpp"

/** 
 * Solves Kepler's equation x = u(x) + e*sin(u(x))
 */

bool Kepler()
{
  int np = MPIComm::world().getNProc();

  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Create a mesh. It will be of type BasisSimplicialMesh, and will
   * be built using a PartitionedLineMesher. */
  MeshType meshType = new BasicSimplicialMeshType();
  const double pi = 4.0*atan(1.0);
  MeshSource mesher = new PartitionedLineMesher(0.0, pi, 200*np, meshType);
  Mesh mesh = mesher.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
      
  /* Create unknown and test functions, discretized using first-order
   * Lagrange interpolants */
  Expr u = new UnknownFunction(new Lagrange(2), "u");
  Expr v = new TestFunction(new Lagrange(2), "v");

  /* Create a discrete space, and discretize the function 1.0 on it */
  DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
  Expr u0 = new DiscreteFunction(discSpace, 1.0, "u0");

  /* Create coordinate function */
  Expr x = new CoordExpr(0);

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad = new GaussianQuadrature(4);

     
  /* Define the weak form */
  Expr ecc = new Sundance::Parameter(0.5, "e");
  Expr eqn = Integral(interior, v*(u - ecc*sin(u) - x), quad);
  Expr bc;

  /* Create a Playa NonlinearOperator object */
  NonlinearProblem nlp(mesh, eqn, bc, v, u, u0, vecType);

#ifdef HAVE_CONFIG_H
  ParameterXMLFileReader reader(searchForFile("SolverParameters/nox.xml"));
#else
  ParameterXMLFileReader reader("nox.xml");
#endif
  ParameterList noxParams = reader.getParameters();

  std::cerr << "solver params = " << noxParams << std::endl;

  NOXSolver solver(noxParams);

  int numEcc = 5;
  double finalEcc = 0.4;
  double maxErr = 0.0;
  for (int r=1; r<=numEcc; r++)
  {
    double e = r*finalEcc/((double) numEcc);
    ecc.setParameterValue(e);
    std::cerr << "--------------------------------------------------------- " << std::endl;
    std::cerr << " solving at eccentricity = " << ecc << std::endl;
    std::cerr << "--------------------------------------------------------- " << std::endl;
    // Solve the nonlinear system
    SolverState<double> status = nlp.solve(solver);
    TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
      runtime_error, "solve failed");
          

    int maxTerms = 400;
    Expr exactSoln = x;
    int nUsed = 1;
    for (int n=1; n<=maxTerms; n++)
    {
      double coeff = 2.0/((double) n) * jn(n, n*e);
      exactSoln = exactSoln + coeff*sin(n*x);
      nUsed++;
      if (fabs(coeff) < 1.0e-8) break;
    }
    std::cerr << "used " << nUsed << " terms in Kapteyn series" << std::endl;
    Expr exactSolnDisc = L2Projector(discSpace, exactSoln).project();
    Expr errExpr = Integral(interior, 
      pow(u0-exactSoln, 2),
      new GaussianQuadrature(8));
    double errorSq = evaluateIntegral(mesh, errExpr);
    std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;
    if (sqrt(errorSq) > maxErr) maxErr = sqrt(errorSq);

    /* Write the field in matlab format */
    FieldWriter w = new MatlabWriter("kepler-e" + Teuchos::toString(e) + ".dat");
    w.addMesh(mesh);
    w.addField("eccentric anomaly", new ExprFieldWrapper(u0[0]));
    w.addField("exact", new ExprFieldWrapper(exactSolnDisc));
    w.write();
  }

  double tol = 1.0e-4;
  return SundanceGlobal::checkTest(maxErr, tol);
}

