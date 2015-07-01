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
#include "SundanceProblemTesting.hpp"


class NitschePoisson2DTest : public LPRectTestBase
{
public:
  /** */
  NitschePoisson2DTest(const Array<int>& n, int k, double C)
    : LPRectTestBase(n), k_(k), C_(C) {}

  /** */
  std::string name() const {return "NitschePoisson2D";}

  /** */
  Expr exactSoln() const 
    {
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      const double pi = 4.0*atan(1.0);
      return cos(pi*x/2.0)*sin(pi*y);    
    }

  /** */
  Array<int> pExpected() const {return tuple<int>(k_+1);}

  /** */
  LinearProblem prob(const Mesh& mesh) const
    {
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr grad = List( dx , dy );

      CellFilter allBdry = domain().north() + domain().south() 
        + domain().east() + domain().west();

      BasisFamily L = new Lagrange( k_ );

      Expr u = new UnknownFunction( L , "u" );
      Expr v = new TestFunction( L , "v" );
      QuadratureFamily quad = new GaussianQuadrature( 2 * k_ );
      
      const double pi = 4.0*atan(1.0);
      Expr force = 1.25*pi*pi*cos(pi*x/2.0)*sin(pi*y);
      Expr eqn = Integral( interior(), 
        (grad*v) * (grad*u) - force * v , quad )
        + NitschePoissonDirichletBC(2, allBdry, quad, 
          1.0, v, u, exactSoln(), C_);
      

      Expr bc;

      return LinearProblem(mesh, eqn, bc, v, u, vecType());
    }

  /**
   * Specify which linear solvers are to be used.
   */
  Array<LPTestSpec> specs() const
    {
      double tol = 0.05;
      return tuple(
        LPTestSpec("amesos.xml", tol, makeSet<int>(1))/*,
                                                        LPTestSpec("aztec-ifpack.xml", tol), 
                                                        LPTestSpec("aztec-ml.xml", tol),
                                          LPTestSpec("belos-ml.xml", tol),
                                          LPTestSpec("bicgstab.xml", tol) */
        );
    }

  void postRunCallback(int meshID, const Mesh& mesh,
    const string& solverFile,
    const Expr& soln) const
    {
      string solverName = StrUtils::before(solverFile, ".");
      string file = "NitschePoisson2D-mesh-" + Teuchos::toString(meshID)
        + "-" + solverName + "-p-" + Teuchos::toString(k_)
        + "-C-" + Teuchos::toString(C_);

      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType());
      Expr uEx = exactSoln();
      L2Projector proj(discSpace, soln - uEx);
      Expr err = proj.project();

      FieldWriter w = new VTKWriter(file);
      w.addMesh(mesh);
      w.addField("soln", new ExprFieldWrapper(soln));
      w.addField("error", new ExprFieldWrapper(err));
      w.write();
    }
private:
  int k_;
  double C_;
  
};


int main( int argc , char **argv )
{

  try
  {
    Sundance::init(&argc, &argv);
    Tabs::showDepth() = false;
    LinearSolveDriver::solveFailureIsFatal() = false;
    int np = MPIComm::world().getNProc();

    LPTestSuite tests;

    Array<int> nx = tuple(4, 8, 16, 32, 64);
    double C = 100.0;

    for (int p=1; p<=3; p++)
    {
      tests.registerTest(rcp(new NitschePoisson2DTest(nx, p, C)));
    }

    bool pass = tests.run();

    Out::root() << "total test status: " << pass << std::endl;

    Sundance::passFailTest(pass);
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); 
  return Sundance::testStatus(); 
}
