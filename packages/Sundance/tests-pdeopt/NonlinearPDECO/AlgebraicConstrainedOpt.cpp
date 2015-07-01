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
#include "PDEOptNonlinearPDEConstrainedObj.hpp"

#include "PlayaBasicLMBFGS.hpp"
#include "PlayaOptBuilder.hpp"


int main(int argc, char** argv)
{
  try
		{
			Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();
      
      int nx = 2;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
      Mesh mesh = mesher.getMesh();
      CellFilter interior = new MaximalCellFilter();
      
      /* Create a vector space factory, used to 
       * specify the low-level linear algebra representation */
      VectorType<double> vecType = new EpetraVectorType();
  
      /* create a discrete space on the mesh */
      DiscreteSpace discreteSpace(mesh, new Lagrange(1), vecType);

      /* initialize the design, state, and multiplier vectors */
      Expr alpha0 = new DiscreteFunction(discreteSpace, 0.5, "alpha0");
      Expr u0 = new DiscreteFunction(discreteSpace, 0.5, "u0");
      Expr lambda0 = new DiscreteFunction(discreteSpace, 0.5, "lambda0");

      /* create symbolic objects for test and unknown functions */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr lambda = new UnknownFunction(new Lagrange(1), "lambda");
      Expr alpha = new UnknownFunction(new Lagrange(1), "alpha");

      /* create quadrature rules of different orders */
      QuadratureFamily q1 = new GaussianQuadrature(1);
      QuadratureFamily q2 = new GaussianQuadrature(2);
      QuadratureFamily q4 = new GaussianQuadrature(4);

      /* Regularization weight */
      double R = 1.0;

      /* Exact solution (computed by Mathematica to 6 digits) */
      double U0 = 1.20404;
      double A0 = -0.20863;
      Expr aExact = new DiscreteFunction(discreteSpace, A0, "alpha0");
      Expr uExact = new DiscreteFunction(discreteSpace, U0, "u0");

      /* Form objective function */
      Expr reg = Integral(interior, 0.5 * R * alpha*alpha, q2);
      Expr fit = Integral(interior, 0.5 * pow(u-1.0, 2.0), q2);

      Expr constraintEqn = Integral(interior, 
        lambda*(alpha - sin(u-sqrt(2.0))), q4);
      Expr L = reg + fit + constraintEqn;

      Expr constraintBC;

      Functional Lagrangian(mesh, L, constraintBC, vecType);
      
      LinearSolver<double> adjSolver 
        = LinearSolverBuilder::createSolver("amesos.xml");
      ParameterXMLFileReader reader("nox.xml");
      ParameterList noxParams = reader.getParameters();
      NOXSolver nonlinSolver(noxParams);

      RCP<ObjectiveBase> obj = rcp(new NonlinearPDEConstrainedObj(
        Lagrangian, u, u0, lambda, lambda0, alpha, alpha0,
        nonlinSolver, adjSolver));

      Vector<double> xInit = obj->getInit();

      bool doFDCheck = true;
      if (doFDCheck)
      {
        Out::root() << "Doing FD check of gradient..." << endl;
        bool fdOK = obj->fdCheck(xInit, 1.0e-6, 2);
        if (fdOK) 
        {
          Out::root() << "FD check OK" << endl;
        }
        else
        {
          Out::root() << "FD check FAILED" << endl;
        }
      }

      RCP<UnconstrainedOptimizerBase> opt 
          = OptBuilder::createOptimizer("basicLMBFGS.xml");
      opt->setVerb(3);

      OptState state = opt->run(obj, xInit);

      bool ok = true;
      if (state.status() != Opt_Converged)
      {
        Out::root()<< "optimization failed: " << state.status() << endl;
        ok = false;
      }
      else
      {
        Out::root() << "opt converged: " << state.iter() << " iterations"
                    << endl;
        Out::root() << "exact solution: U0=" << U0 << endl;
      }
      FieldWriter w = new MatlabWriter("AlgebraicConstrained1D.dat");
      w.addMesh(mesh);
      w.addField("u", new ExprFieldWrapper(u0));
      w.addField("alpha", new ExprFieldWrapper(alpha0));
      w.addField("lambda", new ExprFieldWrapper(lambda0));
      w.write();

          
      double uErr = L2Norm(mesh, interior, u0-uExact, q4);
      double aErr = L2Norm(mesh, interior, alpha0-aExact, q4);
      Out::root() << "error in u = " << uErr << endl;
      Out::root() << "error in alpha = " << aErr << endl;

      double tol = 0.001;
      Sundance::passFailTest(uErr + aErr, tol);
      
      
    }
	catch(std::exception& e)
		{
      cerr << "main() caught exception: " << e.what() << endl;
		}
	Sundance::finalize();
  return Sundance::testStatus(); 
}
