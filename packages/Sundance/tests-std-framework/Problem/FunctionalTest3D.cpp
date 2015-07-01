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
#include "SundanceFunctional.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"

using Sundance::List;
/** 
 *
 */

int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Read the mesh */
      int nx = 4;
      int ny = 2;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, np,
                                                         0.0, 1.0, ny, 1,
                                                         meshType);
      Mesh mesh2D = mesher.getMesh();

      MeshTransformation extruder = new ExtrusionMeshTransformation(0.0, 1.0, 2, meshType);

      Mesh mesh = extruder.apply(mesh2D);



      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      
      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad4 = new GaussianQuadrature(6);

      /* Compute an integral over a fixed integrand */
      const double pi = 4.0*atan(1.0);

      Expr I1 = Integral(interior, x*sin(pi*x), quad4);
      double f1 = evaluateIntegral(mesh, I1);
      cout << "integral of x sin(pi*x) = " << f1 << std::endl;
      double I1Exact = 1.0/pi;
      cout << "exact: " << I1Exact << std::endl;

      double error = fabs(f1 - I1Exact);
      cout << "error = " << fabs(f1 - I1Exact) << std::endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();

      Expr I2 = Integral(interior, x*x*sin(pi*x), quad4);
      double f2 = evaluateIntegral(mesh, I2);
      cout << "integral of x^2 sin(pi*x) = " << f2 << std::endl;
      double I2Exact = (1.0 - 4.0/pi/pi)/pi;
      cout << "exact: " << I2Exact << std::endl;

      error = max(error, fabs(f2 - I2Exact));
      cout << "error = " << fabs(f2 - I2Exact) << std::endl; 

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();

      Expr I3 = Integral(interior, sin(pi*x)*sin(pi*x), quad4);
      double f3 = evaluateIntegral(mesh, I3);
      cout << "integral of sin^2(pi*x) = " << f3 << std::endl;
      double I3Exact = 0.5;
      cout << "exact: " << I3Exact << std::endl;

      error = max(error, fabs(f3 - I3Exact));
      cout << "error = " << fabs(f3 - I3Exact) << std::endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();

      Expr I4 = Integral(interior, x*x*(pi-x)*(pi-x), quad4);
      double f4 = evaluateIntegral(mesh, I4);
      cout << "integral of x^2 (pi-x)^2 = " << f4 << std::endl;
      double I4Exact = pi*pi/3.0 - pi/2 + 1.0/5.0;
      cout << "exact: " << I4Exact << std::endl;

      error = max(error, fabs(f4 - I4Exact));
      cout << "error = " << fabs(f4 - I4Exact) << std::endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();


      Expr I5 = Integral(interior, 1.0, quad4);
      double f5 = evaluateIntegral(mesh, I5);
      cout << "integral of 1.0 = " << f5 << std::endl;
      double I5Exact = 1.0;
      cout << "exact: " << I5Exact << std::endl;

      error = max(error, fabs(f5 - I5Exact));
      cout << "error = " << fabs(f5 - I5Exact) << std::endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();


      /* now compute a functional at a particular value of a field */
      Expr alpha = new UnknownFunction(new Lagrange(2), "alpha");
      Expr beta = new TestFunction(new Lagrange(2), "beta");

      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);
      L2Projector projector(discSpace, x*(pi-x));
      Expr alpha0 = projector.project();

      cout << "computing Integral(0.5*pow(alpha-sin(pi*x), 2)) at alpha0=x*(pi-x)"
           << std::endl;
      //#define BLAHBLAH 1
#ifdef BLAHBLAH
      verbosity<Evaluator>() = 5;
      verbosity<SparsitySuperset>() = 5;
      verbosity<EvaluatableExpr>() = 5;
      verbosity<Assembler>() = 5;
      verbosity<ElementIntegral>() = 5;
      EvalVector::shadowOps() = true;
#endif   
      Expr g = Integral(interior, 0.5*pow(alpha-sin(pi*x), 2.0) , quad4);
      //    Expr g = Integral(interior, 0.5*(alpha-sin(pi*x))*(alpha-sin(pi*x)) , quad4);
      //Expr g = Integral(interior, sin(alpha), quad4);
      //Expr dg = Integral(interior, alpha*beta+cos(alpha0)*beta, quad4);
      Expr dg = Integral(interior, alpha*beta 
                         + (alpha0-sin(pi*x))*beta , quad4);
      Expr bc;
      LinearProblem phony(mesh, dg, bc, beta, alpha, vecType);
      Vector<double> dgVec = phony.getRHS()[0];

      double gExact = 0.5*(-2.0*pi*I1Exact + I3Exact + 2.0*I2Exact + I4Exact);

      //#define BLAHBLAH 1
#ifdef BLAHBLAH
      verbosity<Evaluator>() = 5;
      verbosity<SparsitySuperset>() = 5;
      verbosity<EvaluatableExpr>() = 5;
      verbosity<Assembler>() = 5;
      verbosity<ElementIntegral>() = 5;
      EvalVector::shadowOps() = true;
#endif   

      Functional G(mesh, g, vecType);

      FunctionalEvaluator gEval = G.evaluator(alpha, alpha0);

      cout << "computing function value: " << std::endl;
      double gVal = gEval.evaluate();
      cout << "integral value = " << gVal << std::endl;
      cout << "exact value = " << gExact << std::endl;
      error = max(error, fabs(gVal - gExact));
      cout << "error = " << fabs(gVal - gExact) << std::endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();
      
      /* now compute the derivative of a functional wrt a field variable */
 
      cout << "computing function value and gradient together: " << std::endl;
      Expr dGdAlpha = gEval.evalGradient(gVal);
      cout << "getting vector " << std::endl;
      Vector<double> dgNumVec 
        = DiscreteFunction::discFunc(dGdAlpha)->getVector();
      Vector<double> dgDiff = dgVec - dgNumVec;
      cout << "grad diff = " << std::endl << dgDiff.norm2() << std::endl; 

      cout << "integral value = " << gVal << std::endl;
      error = max(error, fabs(gVal - gExact));
      cout << "error = " << fabs(gVal - gExact) << std::endl;

      MPIComm::world().synchronize();
      MPIComm::world().synchronize();
      cout << "*********************** FD check ***************************** " << std::endl;
      double h = 1.0e-2;
      double diffErr = gEval.fdGradientCheck(h);

      error = max(error, fabs(diffErr));
      cout << "max error = " << error << std::endl;

      double tol = 1.0e-8;
      Sundance::passFailTest(error, tol);


    }
	catch(std::exception& e)
		{
      cout << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
