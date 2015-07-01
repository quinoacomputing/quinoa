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
 * Solves the Poisson equation in 2D
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})

CELL_PREDICATE(CornerPointTest, {return fabs(x[1]) < 1.0e-10 && fabs(x[0])<1.0e-10;})


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
      int nx = 8;
      int ny = 8;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, np,
                                                         0.0, 1.0, ny, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter bdry = new BoundaryCellFilter();
      CellFilter nodes = new DimensionalCellFilter(0);
      CellFilter corner = nodes.subset(new CornerPointTest());

      
      /* Unknown and test functions, using Taylor-Hood discretization */
      Expr ux = new UnknownFunction(new Lagrange(2), "u_x");
      Expr vx = new TestFunction(new Lagrange(2), "v_x");
      Expr uy = new UnknownFunction(new Lagrange(2), "u_y");
      Expr vy = new TestFunction(new Lagrange(2), "v_y");
      Expr p = new UnknownFunction(new Lagrange(1), "p");
      Expr q = new TestFunction(new Lagrange(1), "q");
      Expr u = List(ux, uy);
      Expr v = List(vx, vy);

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

      double pi = 4.0*atan(1.0);
      Expr sx = sin(pi*x);
      Expr cx = cos(pi*x);
      Expr sy = sin(pi*y);
      Expr cy = cos(pi*y);
      Expr psiExact = pow(pi, -3.0) * sx*sy;
      Expr uExact = pow(pi, -2.0)*List(-sx*cy, cx*sy);
      Expr fy = 4.0*cx*sy;

      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy)  - p*(dx*vx+dy*vy)
                          + q*(dx*ux+dy*uy) - vy*fy,
                          quad2);
        
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(bdry, v*(u-uExact), quad4)
        + EssentialBC(corner, q*p, quad4);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, List(vx, vy, q),
                         List(ux, uy, p), vecType);


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/aztec-native.xml"));
#else
      ParameterXMLFileReader reader("aztec-native.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver);


      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("TaylorHoodStokes2d");
      w.addMesh(mesh);
      w.addField("ux", new ExprFieldWrapper(soln[0]));
      w.addField("uy", new ExprFieldWrapper(soln[1]));
      w.addField("p", new ExprFieldWrapper(soln[2]));
      w.write();

      Expr err = List(soln[0], soln[1]) - uExact;
      Expr errExpr = Integral(interior, 
                              err*err,
                              quad4);

      FunctionalEvaluator errInt(mesh, errExpr);

      double errorSq = errInt.evaluate();
      std::cerr << "velocity error norm = " << sqrt(errorSq) << std::endl << std::endl;

      Sundance::passFailTest(sqrt(errorSq), 1.0e-3);

    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
