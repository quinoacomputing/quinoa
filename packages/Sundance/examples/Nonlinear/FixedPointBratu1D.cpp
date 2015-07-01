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

/* 
 * Solve the Bratu problem in 1D using fixed-point iteration 
 */

int main(int argc, char** argv)
{
  try
  {
    int nx = 32;
    double convTol = 1.0e-8;
    double lambda = 0.5;
    Sundance::setOption("nx", nx, "Number of elements");
    Sundance::setOption("tol", convTol, "Convergence tolerance");
    Sundance::setOption("lambda", lambda, "Lambda (parameter in Bratu's equation)");

    Sundance::init(&argc, &argv);

    Out::root() << "Bratu problem (lambda=" << lambda << ")" << endl;
    Out::root() << "Fixed-point iteration" << endl << endl;

    VectorType<double> vecType = new EpetraVectorType();

    MeshType meshType = new BasicSimplicialMeshType();
    MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
    Mesh mesh = mesher.getMesh();

    CellFilter interior = new MaximalCellFilter();
    CellFilter sides = new DimensionalCellFilter(mesh.spatialDim()-1);
    CellFilter left = sides.subset(new CoordinateValueCellPredicate(0, 0.0));
    CellFilter right = sides.subset(new CoordinateValueCellPredicate(0, 1.0));
    
    BasisFamily basis = new Lagrange(1);
    Expr u = new UnknownFunction(basis, "u");
    Expr v = new TestFunction(basis, "v");

    Expr grad = gradient(1);

    Expr x = new CoordExpr(0);



    const double pi = 4.0*atan(1.0);
    Expr uExact = sin(pi*x);
    Expr R = pi*pi*uExact - lambda*exp(uExact);

    QuadratureFamily quad4 = new GaussianQuadrature(4);
    QuadratureFamily quad2 = new GaussianQuadrature(2);

    DiscreteSpace discSpace(mesh, basis, vecType);
    Expr uPrev = new DiscreteFunction(discSpace, 0.5);
    Expr uCur = copyDiscreteFunction(uPrev);

    Expr eqn 
      = Integral(interior, (grad*u)*(grad*v) - v*lambda*exp(uPrev) - v*R, quad4);

    Expr h = new CellDiameterExpr();
    Expr bc = EssentialBC(left+right, v*u/h, quad4); 

    LinearProblem prob(mesh, eqn, bc, v, u, vecType);

    Expr normSqExpr = Integral(interior, pow(u-uPrev, 2.0), quad2);
    Functional normSqFunc(mesh, normSqExpr, vecType);
    FunctionalEvaluator normSqEval = normSqFunc.evaluator(u, uCur);

    LinearSolver<double> linSolver 
      = LinearSolverBuilder::createSolver("amesos.xml");

    Out::root() << "Fixed-point iteration" << endl;
    int maxIters = 20;
    Expr soln ;
    bool converged = false;

    for (int i=0; i<maxIters; i++)
    {
      /* solve for the next u */
      prob.solve(linSolver, uCur);
      /* evaluate the norm of (uCur-uPrev) using 
       * the FunctionalEvaluator defined above */
      double deltaU = sqrt(normSqEval.evaluate());
      Out::root() << "Iter=" << setw(3) << i << " ||Delta u||=" << setw(20)
                  << deltaU << endl; 
      /* check for convergence */  
      if (deltaU < convTol) 
      {
        soln = uCur;
        converged = true;
        break;
      }
      /* get the vector from the current discrete function */
      Vector<double> uVec = getDiscreteFunctionVector(uCur);
      /* copy the vector into the previous discrete function */ 
      setDiscreteFunctionVector(uPrev, uVec);
    } 
    TEUCHOS_TEST_FOR_EXCEPTION(!converged, std::runtime_error, 
      "Fixed point iteration did not converge after " 
      << maxIters << " iterations");
    
    FieldWriter writer = new DSVWriter("FixedPointBratu.dat");
    writer.addMesh(mesh);
    writer.addField("soln", new ExprFieldWrapper(soln[0]));
    writer.write();

    Out::root() << "Converged!" << endl << endl;

    double L2Err = L2Norm(mesh, interior, soln-uExact, quad4);
    Out::root() << "L2 Norm of error: " << L2Err << endl;
    
    Sundance::passFailTest(L2Err, 1.5/((double) nx*nx));
  }
	catch(std::exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
}

