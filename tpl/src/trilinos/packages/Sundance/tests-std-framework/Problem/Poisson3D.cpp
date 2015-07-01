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
#include "SundanceLabelCellPredicate.hpp"
#include "SundanceExodusNetCDFMeshReader.hpp"
#include "SundanceExodusMeshReader.hpp"
#include "SundanceElementIntegral.hpp"

using Sundance::List;
/** 
 * Solves the Poisson equation in 3D
 */



class PolyFunc : public PointwiseUserDefFunctor0
{
public:
  PolyFunc(int n) : PointwiseUserDefFunctor0("P_" + Teuchos::toString(n), 1, 1), n_(n){}

  /** */
  void eval0(const double* vars, double* f) const ;

private:
  int n_;
};


Expr Poly(int n, const Expr& x)
{
  return  new UserDefOp(x, rcp(new PolyFunc(n)));
}



void PolyFunc::eval0(const double* vars, double* f) const
{
  double y = 1.0;
  double x = vars[0];
  for (int i=0; i<n_; i++)
  {
    double t = 1.0;
    for (int j=0; j<n_; j++) t = t*x;
    y = y + 2.0*t - t - t;
  }
  f[0] = y;
}




#if defined(HAVE_SUNDANCE_EXODUS) && defined(Trilinos_DATA_DIR)

int main(int argc, char** argv)
{
  try
		{
      int depth = 0;
      bool useCCode = false;
      Sundance::ElementIntegral::alwaysUseCofacets() = true;
      Sundance::clp().setOption("depth", &depth, "expression depth");
      Sundance::clp().setOption("C", "symb", &useCCode, "Code type (C or symbolic)");
      Sundance::init(&argc, &argv);

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Read the mesh */
      MeshType meshType = new BasicSimplicialMeshType();

      MeshSource mesher 
        = new ExodusMeshReader("cube-0.1", meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter faces = new DimensionalCellFilter(2);
      CellFilter side1 = faces.labeledSubset(1);
      CellFilter side2 = faces.labeledSubset(2);
      CellFilter side3 = faces.labeledSubset(3);
      CellFilter side4 = faces.labeledSubset(4);
      CellFilter side5 = faces.labeledSubset(5);
      CellFilter side6 = faces.labeledSubset(6);

      
      /* Create unknown and test functions, discretized using second-order
       * Lagrange interpolants */
      Expr u = new UnknownFunction(new Lagrange(1), "u");
      Expr v = new TestFunction(new Lagrange(1), "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);
      Expr grad = List(dx, dy, dz);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form */
      //Expr eqn = Integral(interior, (grad*v)*(grad*u) + v, quad);
      
      Expr coeff = 1.0;
#ifdef FOR_TIMING
      if (useCCode)
      {
        coeff = Poly(depth, x);
      }
      else
      {
        for (int i=0; i<depth; i++)
        {
          Expr t = 1.0;
          for (int j=0; j<depth; j++) t = t*x;
          coeff = coeff + 2.0*t - t - t;
        }
      }
#endif
      Expr eqn = Integral(interior, coeff*(grad*v)*(grad*u) /*+ 2.0*v*/, quad2);

      /* Define the Dirichlet BC */
      Expr exactSoln = x;//(x + 1.0)*x - 1.0/4.0;
      Expr h = new CellDiameterExpr();

      WatchFlag watchBC("watch BCs");
      watchBC.setParam("integration setup", 6);
      watchBC.setParam("integration", 6);
      watchBC.setParam("fill", 6);
      watchBC.setParam("evaluation", 6);
      watchBC.deactivate();

      Expr bc = EssentialBC(side4, v*(u-exactSoln), quad4)
        + EssentialBC(side6, v*(u-exactSoln), quad4, watchBC);

      /* We can now set up the linear problem! */
      LinearProblem prob(mesh, eqn, bc, v, u, vecType);

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(searchForFile("SolverParameters/aztec-ml.xml"));
#else
      ParameterXMLFileReader reader("aztec-ml.xml");
#endif
      ParameterList solverParams = reader.getParameters();
      std::cerr << "params = " << solverParams << std::endl;


      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      Expr soln = prob.solve(solver);

#ifndef FOR_TIMING

      DiscreteSpace discSpace(mesh, new Lagrange(1), vecType);
      L2Projector proj1(discSpace, exactSoln);
      L2Projector proj2(discSpace, soln-exactSoln);
      L2Projector proj3(discSpace, pow(soln-exactSoln, 2.0));
      Expr exactDisc = proj1.project();
      Expr errorDisc = proj2.project();
//      Expr errorSqDisc = proj3.project();

      std::cerr << "writing fields" << std::endl;
      /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Poisson3d");
      w.addMesh(mesh);
      w.addField("soln", new ExprFieldWrapper(soln[0]));
      w.addField("exact soln", new ExprFieldWrapper(exactDisc));
      w.addField("error", new ExprFieldWrapper(errorDisc));
//      w.addField("errorSq", new ExprFieldWrapper(errorSqDisc));
      w.write();

      std::cerr << "computing error" << std::endl;

      Expr errExpr = Integral(interior, 
                              pow(soln-exactSoln, 2.0),
                              new GaussianQuadrature(4));

      double errorSq = evaluateIntegral(mesh, errExpr);
      std::cerr << "error norm = " << sqrt(errorSq) << std::endl << std::endl;
#else
      double errorSq = 1.0;
#endif
      double tol = 1.0e-10;
      Sundance::passFailTest(sqrt(errorSq), tol);
    }
	catch(std::exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize(); return Sundance::testStatus(); 
}


#else


int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy Poisson3D PASSED. Enable exodus to run the actual test" << std::endl;
  Sundance::finalize();
  return 0;
}


#endif
