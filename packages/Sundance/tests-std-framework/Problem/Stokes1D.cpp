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
 * Solves the Stokes equation in 2D
 */


CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})

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
      int nx = 60;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx*np,
                                                    meshType);

      



      Mesh mesh = mesher.getMesh();


      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);

      CellFilter left = points.subset(new LeftPointTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr ux = new UnknownFunction(new Lagrange(1), "u_x");
      Expr vx = new TestFunction(new Lagrange(1), "v_x");
      Expr p = new UnknownFunction(new Lagrange(1), "p");
      Expr q = new TestFunction(new Lagrange(1), "q");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr grad = List(dx);
      Expr x = new CoordExpr(0);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      double beta = 0.02;
      Expr h = new CellDiameterExpr();
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          - p*(dx*vx)
                          - (h*h*beta)*(grad*q)*(grad*p) - q*(dx*ux),
                          quad2);
        
      /* Define the Dirichlet BC */
      Expr uInflow = cos(x);
      Expr bc =  EssentialBC(left, vx*(ux-uInflow) , quad4);


      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, q), 
                         List(ux, p), vecType);


     

      /* Create an Aztec solver */
      std::map<int,int> azOptions;
      std::map<int,double> azParams;

      azOptions[AZ_solver] = AZ_gmres;
      azOptions[AZ_precond] = AZ_Neumann;
      azOptions[AZ_poly_ord] = 4;
      //      azOptions[AZ_subdomain_solve] = AZ_ilu;
      //      azOptions[AZ_graph_fill] = 1;
      azOptions[AZ_max_iter] = 5000;
      azParams[AZ_tol] = 1.0e-12;

      LinearSolver<double> solver = new AztecSolver(azOptions,azParams);

      Expr soln = prob.solve(solver);

      /* Write the field in VTK format */
      FieldWriter w = new MatlabWriter("Stokes1d.dat");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("p", new ExprFieldWrapper(soln[1]));
      w.write();

      Expr uxErr = soln[0] - 1.0;
      Expr errExpr = Integral(interior, 
                              uxErr*uxErr,
                              new GaussianQuadrature(4));

      FunctionalEvaluator errInt(mesh, errExpr);

      double errorSq = errInt.evaluate();
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

      
      double tol = 1.0e-10;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  TimeMonitor::summarize();
  
}
