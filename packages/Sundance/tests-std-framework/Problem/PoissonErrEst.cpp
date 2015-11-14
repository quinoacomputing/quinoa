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
#include "SundanceBubble.hpp"
#include "SundanceRivaraDriver.hpp"


/** 
 * Solves the Poisson equation with error control
 */

class MyProb
{
public:
  CellFilter interior_;
  CellFilter sides_;
  CellFilter outerBdry_;
  Expr uEx_;
  Expr rho_;


public:
  /** constructor */
  MyProb()
  {
    interior_ = new MaximalCellFilter();
    sides_ = new DimensionalCellFilter(2);
    outerBdry_ = sides_.labeledSubset(1);
    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);
    Expr z = new CoordExpr(2);

    Expr r2 = x*x + y*y;
    double a = 100.0;
    Expr b = 1.0 + a*r2;
    uEx_ = log(b/(1+a));
    rho_ = 4.0*a/b/b;
  }
  
 
  
  Expr resid(const Expr& v, const Expr& u)
  {
    /* Create differential operators, gradient, and coordinate functions */
    Expr grad = gradient(3);
  
    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad4 = new GaussianQuadrature(4);
  
    /* Define the weak form of the residual */
    return Integral(interior_, (grad*v)*(grad*u) + v*rho_, quad4);
  }
};



int main(int argc, char** argv)
{
  
  try
    {
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh object to fill in from a file. It will be of type BasisSimplicialMesh. */
      MeshType meshType = new BasicSimplicialMeshType();

      /* read the mesh from the file disk*/
      MeshSource meshReader = new ExodusMeshReader("disk", meshType);
      Mesh mesh = meshReader.getMesh();

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      Expr e = new UnknownFunction(new Bubble(3), "e");
      Expr eHat = new TestFunction(new Bubble(3), "eHat");

      MyProb myProb;

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad4 = new GaussianQuadrature(4);
     
      /* Define the weak form */
      Expr eqn = myProb.resid(v, u);
      /* Define the Dirichlet BC */
      Expr bc = EssentialBC(myProb.outerBdry_, v*u, quad4);

      Expr uEx = myProb.uEx_;

      LinearSolver<double> solver 
	= LinearSolverBuilder::createSolver("amesos.xml");

      int maxTries = 4;
      double goal = 1.0e-3;
      double errNorm = 1.0;
      for (int i=0; i<maxTries; i++)
	{
	  LinearProblem prob(mesh, eqn, bc, v, u, vecType);
	  Expr soln = prob.solve(solver);

	  DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
	  L2Projector errProj(discSpace, soln-uEx);
	  Expr err0 = errProj.project();
	  L2Projector uExProj(discSpace, uEx);
	  Expr u0 = uExProj.project();
	  
	  /* weak form for error estimation */
	  Expr grad = gradient(3);
	  Expr errEqn = Integral(myProb.interior_, (grad*eHat)*(grad*e), quad4) - myProb.resid(eHat, soln);
	  Expr errBC;
	  LinearProblem errProb(mesh, errEqn, errBC, eHat, e, vecType);
	  Expr errEst = errProb.solve(solver);
      
      
	  FieldWriter w = new VTKWriter("PoissonErrEst-" + Teuchos::toString(i));
	  w.addMesh(mesh);
	  w.addField("numerical soln", new ExprFieldWrapper(soln[0]));
	  w.addField("exact soln", new ExprFieldWrapper(u0));
	  w.addField("error", new ExprFieldWrapper(err0));
	  w.addField("error estimate", new ExprFieldWrapper(errEst));
	  w.write();

	  errNorm = L2Norm(mesh, myProb.interior_, soln-uEx, quad4);
	  Out::root() << "error norm = " << errNorm << endl;
	  Out::root() << "goal = " << goal << endl;
	  if (errNorm < goal)
	    {
	      Out::root() << "error norm is within goal!" << endl;
	      break;
	    }

	  RefinementTransformation refiner(meshType, errEst, goal, 1.0e-4);
	  mesh = refiner.apply(mesh);
	}
      double tol = 1.0e-4;
      Sundance::passFailTest(errNorm, tol);
      
    }
  catch(std::exception& e)
    {
      std::cerr << e.what() << std::endl;
    }
  Sundance::finalize(); return Sundance::testStatus(); 
}


