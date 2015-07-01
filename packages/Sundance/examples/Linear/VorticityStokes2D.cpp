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
#include "Teuchos_XMLParameterListWriter.hpp"

using Sundance::List;
/** 
 * Solves the 2D Stokes driven cavity in the vorticity-streamfunction 
 * formulation
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})


int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      int nx = 32;
      int ny = 32;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, np,
                                                         0.0, 1.0, ny, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter left = edges.coordSubset(0,0.0);
      CellFilter right = edges.coordSubset(0,1.0);
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

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      Expr eqn = Integral(interior, (grad*vPsi)*(grad*psi) 
                          + (grad*vOmega)*(grad*omega) + vPsi*omega, quad2)
        + Integral(top, -1.0*vPsi, quad4);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bottom, vOmega*psi, quad2) 
        + EssentialBC(top, vOmega*psi, quad2) 
        + EssentialBC(left, vOmega*psi, quad2) 
        + EssentialBC(right, vOmega*psi, quad2);



      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vPsi, vOmega), 
                         List(psi, omega), vecType);


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver("bicgstab.xml");


      std::cerr << "starting solve..." << std::endl;

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("VorticityStokes2D");
      w.addMesh(mesh);
      w.addField("psi", new ExprFieldWrapper(soln[0]));
      w.addField("omega", new ExprFieldWrapper(soln[1]));
      w.write();

      /* As a check, we integrate the vorticity over the domain. By 
       * Stokes' theorem this should be equal to the line integral
       * of the velocity around the boundary (which is 1). */
      Expr totalVorticityExpr = Integral(interior, soln[1], quad2);
      double totalVorticity = evaluateIntegral(mesh, totalVorticityExpr);
      std::cerr << "total vorticity = " << totalVorticity << std::endl;

      double tol = 1.0e-4;
      Sundance::passFailTest(fabs(totalVorticity-1.0), tol);
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
