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
#include "SundancePointwiseUserDefFunctor.hpp"
#include "SundanceUserDefOp.hpp"
#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_Playa_Group.hpp"


/** Evaluate u*v - 6.0 */
class F1 : public PointwiseUserDefFunctor1
{
public:
  F1() : PointwiseUserDefFunctor1("F1", 2, 1){;}
  virtual ~F1(){;}
  void eval1(const double* vars, double* f, double* df) const ;
  void eval0(const double* vars, double* f) const ;
};


void F1::eval1(const double* vars, double* f, double* df) const
{
  f[0] = vars[0]*vars[1] - 6.0;
  df[0] = vars[1];
  df[1] = vars[0];
}

void F1::eval0(const double* vars, double* f) const
{
  f[0] = vars[0]*vars[1] - 6.0;
}


/** Evaluate u^2 - v - 1.0 */
class F2 : public PointwiseUserDefFunctor1
{
public:
  F2() : PointwiseUserDefFunctor1("F2", 2, 1){;}
  virtual ~F2(){;}
  void eval1(const double* vars, double* f, double* df) const ;
  void eval0(const double* vars, double* f) const ;
};


void F2::eval1(const double* vars, double* f, double* df) const
{
  f[0] = vars[0]*vars[0] - vars[1] - 1.0;
  df[0] = 2.0*vars[0];
  df[1] = -1.0;
}

void F2::eval0(const double* vars, double* f) const
{
  f[0] = vars[0]*vars[0] - vars[1] - 1.0;
}






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
      MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 1*np, meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();

      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      Expr u1 = new UnknownFunction(new Lagrange(1), "u1");
      Expr v1 = new TestFunction(new Lagrange(1), "v1");
      Expr u2 = new UnknownFunction(new Lagrange(1), "u2");
      Expr v2 = new TestFunction(new Lagrange(1), "v2");

      /* Create a discrete space, and discretize the function 1.5 on it */
      BasisFamily L1 = new Lagrange(1);
      DiscreteSpace discSpace(mesh, Sundance::List(L1, L1), vecType);
      Expr u0 = new DiscreteFunction(discSpace, 1.5, "u0");

      Expr f1 = new UserDefOp(List(u1, u2), rcp(new F1()));
      Expr f2 = new UserDefOp(List(u1, u2), rcp(new F2()));

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad = new GaussianQuadrature(2);

      /* Now we set up the weak form of our equation. */
      Expr eqn = Integral(interior, v1*f1 + v2*f2, quad);

      /* There are no boundary conditions for this problem, so the
       * BC expression is empty */
      Expr bc;

        /* We can now set up the nonlinear problem! */
        NonlinearProblem prob(mesh, eqn, bc, List(v1, v2), List(u1, u2), u0, vecType);


#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/nox.xml"));
#else
      ParameterXMLFileReader reader("nox.xml");
#endif
      ParameterList noxParams = reader.getParameters();

      std::cerr << "solver params = " << noxParams << std::endl;

      NOXSolver solver(noxParams);

      prob.solve(solver);

      Expr errExpr = Integral(interior, 
                              pow(u0[0]-2.0, 2) + pow(u0[1]-3.0, 2),
                              new GaussianQuadrature(2));

      double errorSq = evaluateIntegral(mesh, errExpr);
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;

      double tol = 1.0e-8;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(std::exception& e)
		{
      std::cerr << e.what() << std::endl;
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}
