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
 * Solves the heat equation in 1D using Crank-Nicolson timestepping
 */

const double pi = 4.0*atan(1.0);

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-pi) < 1.0e-10;})

int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 512;
      MeshSource mesher = new PartitionedLineMesher(0.0, pi, nx*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter points = new DimensionalCellFilter(0);
      CellFilter leftPoint = points.subset(new LeftPointTest());
      CellFilter rightPoint = points.subset(new RightPointTest());

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      BasisFamily bas = new Lagrange(2);
      Expr u = new UnknownFunction(bas, "u");
      Expr v = new TestFunction(bas, "v");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr x = new CoordExpr(0);

      /* The initial profile is u(x,0)=sin(x). Project this onto a discrete function */
      DiscreteSpace discSpace(mesh, bas, vecType);
      L2Projector projector(discSpace, sin(x));
      Expr u0 = projector.project();

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2);

      int nSteps = 20;
      double dt = 0.2/((double) nSteps);


      /* Define the weak form, semidiscretized in time */
      Expr eqn = Integral(interior, v*(u-u0) + dt/2.0*(dx*v)*((dx*u)+(dx*u0)), quad); 
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(leftPoint, v*u, quad) + EssentialBC(rightPoint, v*u,quad);

      /* We can now set up the linear problem! */
      LinearProblem stepProb(mesh, eqn, bc, v, u, vecType); 

      ParameterXMLFileReader reader("aztec-ifpack.xml");
      ParameterList solverParams = reader.getParameters();

      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      /* define an expression for the exact solution */
      Expr t = new Sundance::Parameter(0.0);
      Expr uExact = exp(-t) * sin(x);
      Expr errExpr = Integral(interior, 
                              pow(u0-uExact, 2),
                              new GaussianQuadrature(6));

      /* loop over timesteps */
      for (int i=0; i<nSteps; i++)
      {
        double tVal = (i+1)*dt;
        t.setParameterValue(tVal);
        Expr uNext = stepProb.solve(solver);
        updateDiscreteFunction(uNext, u0);
        double errorSq = evaluateIntegral(mesh, errExpr);
        FieldWriter writer = new MatlabWriter("cnheat1D-" + Teuchos::toString(i+1));
        writer.addMesh(mesh);
        writer.addField("u", new ExprFieldWrapper(u0[0]));
        writer.write();
        Out::os() << "t=" << tVal << " error norm = " << sqrt(errorSq) << std::endl;
      }
      double errorSq = evaluateIntegral(mesh, errExpr);
      double tol = 1.0e-4;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
