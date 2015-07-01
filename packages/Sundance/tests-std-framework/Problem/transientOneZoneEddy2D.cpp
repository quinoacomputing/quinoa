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
#include "SundanceZeroExpr.hpp"



Expr rot(const Expr& s)
{
  Expr nabla = gradient(2);
  return List(nabla[1]*s, -nabla[0]*s);
}

Expr sqr(const Expr& x) {return x*x;}



/** 
 * Solves the heat equation in 1D using Crank-Nicolson timestepping
 */

const double pi = 4.0*atan(1.0);

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-pi) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;}) 
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-pi) < 1.0e-10;}) 



int main(int argc, char** argv)
{
  try
		{
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      int nx = 32;
      int ny = 32;
      MeshSource mesher = new PartitionedRectangleMesher(0.0, pi, nx, 1, 
        0.0,  pi, ny, 1, meshType);
      Mesh mesh = mesher.getMesh();


      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());
      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());


      Expr mu1 = 1.0;
      Expr sigma1 = 1.0;
      Expr h = new CellDiameterExpr();

      /* Create unknown function */
      Expr H1x = new UnknownFunction(new Lagrange(1), "H1x");
      Expr H1y = new UnknownFunction(new Lagrange(1), "H1y");
      Expr H1Next = List(H1x, H1y);

      Expr E1Next = new UnknownFunction(new Lagrange(1), "E1z");

      /* Create differential operator and coordinate function */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* Project the initial value onto a discrete function */
      Expr E0 = sin(x)*sin(y);
      BasisFamily L1 = new Lagrange(1);
      DiscreteSpace discSpace(mesh, List(L1, L1, L1), vecType);
      L2Projector projector(discSpace, List(dy*E0, -dx*E0, E0));
      Expr u0 = projector.project();
      Expr H1Prev = List(u0[0], u0[1]);
      Expr E1Prev = u0[2];

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      int nSteps = 100;
      double dt = 0.1/((double) nSteps);
      Expr t = new Sundance::Parameter(0.0);
      Expr tPrev = new Sundance::Parameter(0.0);

      Expr HBdryNext = List(sin(x)*cos(y), -cos(x)*sin(y))*exp(-t);
      Expr HBdryPrev = List(sin(x)*cos(y), -cos(x)*sin(y))*exp(-tPrev);
      Expr HBdryAvg = (HBdryNext + HBdryPrev)/2.0;
      
      Expr H1Avg = (H1Next + H1Prev)/2.0;
      Expr E1Avg = (E1Next + E1Prev)/2.0;
      Expr B1Next = mu1*H1Next;
      Expr B1Avg = mu1*H1Avg;
      Expr B1Dot = mu1*(H1Next - H1Prev)/dt;

      Expr sqResid = 
        Integral(interior, 
          sqr(curl(H1Avg)-sigma1*E1Avg) + sqr(rot(E1Avg)+B1Dot) 
          + sqr(div(B1Avg)), quad)
        + Integral(left+right+top+bottom, sqr(HBdryAvg-H1Avg), quad4);
      
      Functional sq(mesh, sqResid, vecType);
      
      Expr dum;
      Expr zero = new Sundance::ZeroExpr();
      Expr u = List(H1x, H1y, E1Next);
      Expr v0 = List(zero,zero,zero);
      LinearProblem stepProb = sq.linearVariationalProb(u, v0, u, dum, dum);


      ParameterXMLFileReader reader("amesos.xml");
      ParameterList solverParams = reader.getParameters();

      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr errExpr = Integral(interior, 
        sqr(u0[2]-E0*exp(-t)), new GaussianQuadrature(6));

      double errorSq = evaluateIntegral(mesh, errExpr);
      Out::os() << "initial error = " << sqrt(errorSq) << std::endl;

      /* loop over timesteps */
      for (int i=1; i<=nSteps; i++)
      {

        double tPrevVal = (i-1)*dt;
        double tVal = i*dt;
        tPrev.setParameterValue(tPrevVal);
        t.setParameterValue(tVal);
        
        Expr uNext = stepProb.solve(solver);
        updateDiscreteFunction(uNext, u0);

        double errorSq = evaluateIntegral(mesh, errExpr);

        Out::os() << "step=" << i << " error=" << sqrt(errorSq) << std::endl;

        FieldWriter writer = new VTKWriter("eddy2D-" + Teuchos::toString(i));
        writer.addMesh(mesh);
        writer.addField("Hx", new ExprFieldWrapper(u0[0]));
        writer.addField("Hy", new ExprFieldWrapper(u0[1]));
        writer.addField("Ez", new ExprFieldWrapper(u0[2]));
        writer.write();
      }
      errorSq = evaluateIntegral(mesh, errExpr);
      double tol = 1.0e-4;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
