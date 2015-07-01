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

#include "PlayaNOXSolver.hpp"

using Sundance::List;
/** 
 * Solves the Navier-Stokes equations on the lid-driver cavity
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int n=16;
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, n, np,
                                                         0.0, 1.0, n, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr psi = new UnknownFunction(new Lagrange(1), "psi");
      Expr vPsi = new TestFunction(new Lagrange(1), "vPsi");
      Expr omega = new UnknownFunction(new Lagrange(1), "omega");
      Expr vOmega = new TestFunction(new Lagrange(1), "vOmega");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* A parameter expression for the Reynolds number */
      Expr reynolds = new Sundance::Parameter(20.0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad1 = new GaussianQuadrature(1);
      QuadratureFamily quad2 = new GaussianQuadrature(2);

      /* Define the weak form */
      Expr psiEqn = Integral(interior, (grad*vPsi)*(grad*psi) + vPsi*omega, 
                             quad2)
        + Integral(top, -1.0*vPsi, quad1);

      Expr omegaEqn = Integral(interior, (grad*vOmega)*(grad*omega)
                               + reynolds*vOmega*((dy*psi)*dx*omega 
                                                  - (dx*psi)*dy*omega),
                               quad1);
      
      Expr eqn = omegaEqn + psiEqn;

      /* Define the Dirichlet BC */
      CellFilter walls = left + bottom + right;
      Expr bc = EssentialBC(walls, vOmega*psi, quad2) 
        + EssentialBC(top, vOmega*psi, quad2) ;

      BasisFamily L1 = new Lagrange(1);
      DiscreteSpace discSpace(mesh, Sundance::List(L1, L1), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 0.0, "u0");

      /* Create a Playa NonlinearOperator object */
      NonlinearProblem prob(mesh, eqn, bc, List(vPsi, vOmega),
        List(psi, omega), u0, vecType);


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/nox.xml"));
#else
      ParameterXMLFileReader reader("nox.xml");
#endif
      ParameterList noxParams = reader.getParameters();

      std::cerr << "solver params = " << noxParams << std::endl;

      NOXSolver solver(noxParams);

      int numReynolds = 10;
      double finalReynolds = 100.0;
      for (int r=0; r<numReynolds; r++)
        {
          double Re = r*finalReynolds/((double) numReynolds-1);
          reynolds.setParameterValue(Re);
          std::cerr << "--------------------------------------------------------- " << std::endl;
          std::cerr << " solving for Reynolds Number = " << Re << std::endl;
          std::cerr << " reynolds = " << reynolds << std::endl;
          std::cerr << "--------------------------------------------------------- " << std::endl;
          // Solve the nonlinear system
          SolverState<double> status = prob.solve(solver);
          TEUCHOS_TEST_FOR_EXCEPTION(status.finalState() != SolveConverged,
            runtime_error, "solve failed");

          /* Write the field in VTK format */
          FieldWriter w = new VTKWriter("vns-r" + Teuchos::toString(Re));
          w.addMesh(mesh);
          w.addField("streamfunction", new ExprFieldWrapper(u0[0]));
          w.addField("vorticity", new ExprFieldWrapper(u0[1]));
          w.write();
        }

      /* As a check, we integrate the vorticity over the domain. By 
       * Stokes' theorem this should be equal to the line integral
       * of the velocity around the boundary. */
      Expr totalVorticityExpr = Integral(interior, u0[1], quad2);
      double totalVorticity = evaluateIntegral(mesh, totalVorticityExpr);
      std::cerr << "total vorticity = " << totalVorticity << std::endl;

      double tol = 1.0e-3;
      Sundance::passFailTest(fabs(totalVorticity-1.0), tol);
      TimeMonitor::summarize();
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}

  
}
