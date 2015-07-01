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
#include "SundanceCubicHermite.hpp"
#include "SundanceProblemTesting.hpp"


/** 
 * This class sets up a test that checks convergence rate for solution
 * of the Poisson equation in 1D with Lagrange(2) basis functions.
 * The test problem is 
 * \f[ u'' = -\frac{1}{4}\pi^2 \sin(\pi x/2) \f]
 * with boundary conditions
 * \f[ u(0)=1, \;\;\; u'(1)=0. \f]
 * The solution is \f$u(x)=\sin(\pi x/2)\f$.
 *
 * 
 */
class SimplePoisson1DTest : public LP1DTestBase
{
public:
  /** 
   * Construct the test object. The constructor calls the base class
   * constructor with a list of mesh sizes to be used. 
   */
  SimplePoisson1DTest(const Array<int>& nx) 
    : LP1DTestBase(nx) {}

  /** Return a descriptive name for the test */
  std::string name() const {return "SimplePoisson1D";}

  /** Return the expected order of accuracy for the solutions. Because a
   * problem might have multiple field variables, the expected orders of 
   * accuracy for the different variables are returned as an array.  */
  Array<int> pExpected() const {return tuple<int>(3);}

  /** 
   * Return an expression for the exact solution to the test problem.
   */
  Expr exactSoln() const 
    {
      Expr x = coord(0);
      const double pi = 4.0*atan(1.0);

      return sin(pi*x/2.0);
    }

  /** 
   * Return a LP for the specified mesh. This function implements the
   * prob() pure virtual function of class LPTestBase.
   */
  LinearProblem prob(const Mesh& mesh) const
    {
      const double pi = 4.0*atan(1.0);
      CellFilter left = domain().left();

      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");
      Expr dx = gradient(1);
      Expr x = coord(0);

      QuadratureFamily quad = new GaussianQuadrature(4);

      Expr eqn = Integral(interior(), (dx*v)*(dx*u), quad)
        + Integral(interior(), -0.25*pi*pi*sin(pi*x/2.0)*v, quad);

      Expr bc = EssentialBC(left, v*u, quad);
      
      return LinearProblem(mesh, eqn, bc, v, u, vecType());
    }

  /**
   * Specify which linear solvers are to be used.
   */
  Array<LPTestSpec> specs() const
    {
      /* Use the standard solvers minus belos-ifpack, which has convergence
       * troubles on the finest mesh */
      double tol = 0.05;
      return tuple(
        LPTestSpec("amesos.xml", tol, makeSet<int>(1)),
        LPTestSpec("aztec-ifpack.xml", tol),
        LPTestSpec("aztec-ml.xml", tol),
        LPTestSpec("belos-ml.xml", tol),
        LPTestSpec("bicgstab.xml", tol)
        );
    }
};




/** Solves the Poisson equation in 1D with user-specified basis functions */
class Poisson1DTest : public LP1DTestBase
{
public:
  /** */
  Poisson1DTest(const Array<int>& nx, const BasisFamily& basis) 
    : LP1DTestBase(nx), basis_(basis) {}

  /** */
  std::string name() const {
    TeuchosOStringStream ss;
    ss << "Poisson1D(basis=";
    basis_.print(ss);
    ss << ")";
    return ss.str();
  }

  /** */
  Array<int> pExpected() const {return tuple<int>(basis_.order()+1);}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);
      const double pi = 4.0*atan(1.0);

      return sin(3.0*pi*x/2.0);
    }

  /** */
  LinearProblem prob(const Mesh& mesh) const
    {
      CellFilter left = domain().left();

      Expr u = new UnknownFunction(basis_, "u");
      Expr v = new TestFunction(basis_, "v");
      Expr dx = gradient(1);
      Expr x = coord(0);

      QuadratureFamily quad = new GaussianQuadrature(2*basis_.order());

      const double pi = 4.0*atan(1.0);
      Expr eqn = Integral(interior(), (dx*v)*(dx*u), quad)
        + Integral(interior(), -2.25*pi*pi*sin(3.0*pi*x/2.0)*v, quad);

      Expr bc = EssentialBC(left, v*u, quad);
      
      
      return LinearProblem(mesh, eqn, bc, v, u, vecType());
    }

  /** */
  Array<LPTestSpec> specs() const
    {
      double tol = 0.05;
      return tuple(
        LPTestSpec("amesos.xml", tol, makeSet<int>(1)),
        LPTestSpec("aztec-ifpack.xml", tol),
        LPTestSpec("aztec-ml.xml", tol),
        LPTestSpec("belos-ml.xml", tol),
        LPTestSpec("bicgstab.xml", tol)
        );
    }
private:
  BasisFamily basis_;
};




/** Performs L2 projection in 1D with user-specified basis functions */
class Projection1DTest : public LP1DTestBase
{
public:
  /** */
  Projection1DTest(const Array<int>& nx, const BasisFamily& basis) 
    : LP1DTestBase(nx), basis_(basis) {}

  /** */
  std::string name() const 
    {
      TeuchosOStringStream ss;
      ss << "Projection1D(basis=";
      basis_.print(ss);
      ss << ")";
      return ss.str();
    }

  /** */
  Array<int> pExpected() const {return tuple<int>(basis_.order()+1);}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);
      return 10.0*x*exp(-1.0*x);
    }

  /** */
  LinearProblem prob(const Mesh& mesh) const
    {
      Expr u = new UnknownFunction(basis_, "u");
      Expr v = new TestFunction(basis_, "v");
      Expr x = coord(0);

      int p = basis_.order();
      QuadratureFamily quad = new GaussianQuadrature(2*p);

      Expr ex = exactSoln();
      Expr eqn = Integral(interior(), v*(u-ex), quad);
      Expr bc;
      
      return LinearProblem(mesh, eqn, bc, v, u, vecType());
    }

  /** */
  Array<LPTestSpec> specs() const
    {
      double tol = 0.05;
      return tuple(
        LPTestSpec("amesos.xml", tol, makeSet<int>(1)),
        LPTestSpec("aztec-ifpack.xml", tol),
        LPTestSpec("aztec-ml.xml", tol),
        LPTestSpec("bicgstab.xml", tol)
        );
    }
private:
  BasisFamily basis_;
};



/** Solves the Helmholtz equation in 1D with Lagrange(2) basis functions */
class Helmholtz1DTest : public LP1DTestBase
{
public:
  /** */
  Helmholtz1DTest(const Array<int>& nx) 
    : LP1DTestBase(0.0, atan(1.0), nx) {}

  /** */
  std::string name() const {return "Helmholtz1D";}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);
      return cos(x) + sin(x);
    }

  /** */
  Array<int> pExpected() const {return tuple(3);}

  /** */
  LinearProblem prob(const Mesh& mesh) const
    {
      CellFilter left = domain().left();
      CellFilter right = domain().right();

      Expr u = new UnknownFunction(new Lagrange(2), "u");
      Expr v = new TestFunction(new Lagrange(2), "v");
      Expr dx = gradient(1);
      Expr x = coord(0);

      QuadratureFamily quad = new GaussianQuadrature(4);

      Expr eqn = Integral(interior(), (dx*v)*(dx*u) - v*u, quad);
      Expr bc = EssentialBC(left, v*(u-cos(x)), quad);
      
      return LinearProblem(mesh, eqn, bc, v, u, vecType());
    }

  /** */
  Array<LPTestSpec> specs() const
    {
      double tol = 0.05;
      return tuple(
        LPTestSpec("amesos.xml", tol, makeSet<int>(1)),
        LPTestSpec("aztec-ifpack.xml", tol),
        LPTestSpec("aztec-ml.xml", tol),
        LPTestSpec("belos-ml.xml", tol),
        LPTestSpec("bicgstab.xml", tol)
        );
    }
private:
};


/** */
class Coupled1DTest : public LP1DTestBase
{
public:
  /** */
  Coupled1DTest(const Array<int>& nx) : LP1DTestBase(nx) {}

  /** */
  std::string name() const {return "Coupled1D";}

  /** */
  Expr exactSoln() const 
    {
      Expr x = coord(0);

      Expr x2 = x*x;
      Expr x3 = x*x2;
      
      Expr uExact = (1.0/120.0)*x2*x3 - 1.0/36.0 * x3 + 7.0/360.0 * x;
      Expr vExact = 1.0/6.0 * x * (x2 - 1.0);

      return List(vExact, uExact);
    }

  /** */
  Array<int> pExpected() const {return tuple(2, 2);}

  /** */
  LinearProblem prob(const Mesh& mesh) const
    {
      CellFilter left = domain().left();
      CellFilter right = domain().right();


      Expr u = new UnknownFunction(new Lagrange(3), "u");
      Expr v = new UnknownFunction(new Lagrange(1), "v");
      Expr du = new TestFunction(new Lagrange(3), "du");
      Expr dv = new TestFunction(new Lagrange(1), "dv");

      Expr dx = gradient(1);
      Expr x = coord(0);


      QuadratureFamily quad = new GaussianQuadrature(10);

      Expr eqn = Integral(interior(), 
        (dx*du)*(dx*u) + du*v + (dx*dv)*(dx*v) + x*dv, 
        quad);

      Expr bc = EssentialBC(left, du*u + dv*v, quad)
        + EssentialBC(right, du*u + dv*v, quad);


      return LinearProblem(mesh, eqn, bc, 
        List(dv,du), List(v,u), vecType());
      
    }  

  /** */
  Array<LPTestSpec> specs() const
    {
      double tol = 0.05;
      return tuple(
        LPTestSpec("amesos.xml", tol, makeSet<int>(1)),
        LPTestSpec("aztec-ifpack.xml", tol),
        LPTestSpec("belos-ifpack.xml", tol),
        LPTestSpec("bicgstab.xml", tol)
        );
    }
};







int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);
    Tabs::showDepth() = false;
    LinearSolveDriver::solveFailureIsFatal() = false;
    int np = MPIComm::world().getNProc();

    LPTestSuite tests;

    Array<int> nx = tuple(8, 16, 24, 32, 40);

    tests.registerTest(rcp(new SimplePoisson1DTest(nx)));

    for (int p=1; p<=3; p++)
    {
      if (np==1 || np % 4 == 0) 
      {
        tests.registerTest(rcp(new Poisson1DTest(nx, new Lagrange(p))));
        tests.registerTest(rcp(new Poisson1DTest(nx, new Bernstein(p))));
        tests.registerTest(rcp(new Projection1DTest(nx, new Lagrange(p))));
        tests.registerTest(rcp(new Projection1DTest(nx, new Bernstein(p))));
      }
      else
      {
        Out::root() << "skipping 1D tests needing numprocs=1 or 4" << std::endl;
      }
    }

    tests.registerTest(rcp(new Helmholtz1DTest(nx))); 

    
    tests.registerTest(rcp(new Coupled1DTest(nx))); 

    bool pass = tests.run();

    Out::root() << "total test status: " << pass << std::endl;

    Sundance::passFailTest(pass);
  }
	catch(std::exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 
}
