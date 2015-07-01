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
#include "SundanceSpectralPreprocessor.hpp"


/** 
 * Solves the Poisson equation in 1D
 */


string expand(const Expr& e)
{
  Tabs tab;
  Array<Array<Expr> > terms;
  SpectralPreprocessor::expandSpectral(e, terms);
  TeuchosOStringStream os;

  os << std::endl << tab << "Terms: " << std::endl;
  for (int i=0; i<terms.size(); i++)
  {
    Tabs tab1;
    os << tab1 << "term=" << i << std::endl;
    for (int j=0; j<terms[i].size(); j++)
    {
      Tabs tab2;
      os << tab2 << "factor " << j << " = " << terms[i][j] << std::endl;
    }
  }

  os << std::endl << tab << "projected expr: " << std::endl;
  {
    Tabs tab1;
    os << tab1 << SpectralPreprocessor::projectSpectral(terms) << std::endl;
  }

  return os.str();
}


Expr LegendreP(int n, const Expr& x)
{
  TEUCHOS_TEST_FOR_EXCEPT(n<0);

  if (n==0) return 1.0;
  if (n==1) return x;
  return ((2*(n-1)+1)*x*LegendreP(n-1,x) - (n-1)*LegendreP(n-2,x))/((double) n);
  
}

bool SpectralPoisson1D()
{
  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Create a mesh. It will be of type BasisSimplicialMesh, and will
   * be built using a PartitionedLineMesher. */
  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, 64, meshType);
  Mesh mesh = mesher.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
  CellFilter points = new DimensionalCellFilter(0);
  CellFilter leftPoint = points.subset(new CoordinateValueCellPredicate(0,0.0));

  /* Create the Spectral Basis */
  int ndim = 1;
  int order = 4;

  SpectralBasis sbasis = new HermiteSpectralBasis(ndim, order); 

  Out::os() << "created the spectral basis" << std::endl;

  /* Create the Spectral Unknown and test functions */

  Expr u = new UnknownFunction(new Lagrange(2),sbasis, "u");
  Expr v = new TestFunction(new Lagrange(2), sbasis, "v");

  Out::os() <<"Unknown and Test Functions " << std::endl; 
  Out::os() << "u=" << u << std::endl;
  Out::os() << "v=" << v << std::endl;

  /* Create differential operator and coordinate function */
  Expr dx = new Derivative(0);
  Expr x = new CoordExpr(0);

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad = new GaussianQuadrature(8);


  /* Define the Stochastic RHS */

  Array<Expr> Q(sbasis.nterms());
  Q[0] = -211.0/100.0;
  Q[1] = -33.0/25.0;
  Q[2] = -16.0/25.0;
  Q[3] = -3.0/50.0;
  Q[4] = -1.0/200.0;
  for (int i=5; i<sbasis.nterms(); i++) Q[i] = 0.0;

  Array<Expr> K(sbasis.nterms());
  Array<Expr> w(sbasis.nterms());
  double p2 = 1.0;
  double p10 = 1.0;
  for (int n=0; n<sbasis.nterms(); n++)
  {
    if (n <= sbasis.nterms()/2) 
    {
      K[n]=1.0/p2;
      w[n]=x*(2.0-x)/p10;
    }
    else 
    {
      K[n]=0.0;
      w[n]=0.0;
    }
    p2 = 2.0*p2;
    p10 = 10.0*p10;
  }

  for (int k=0; k<sbasis.nterms(); k++)
  {
    Tabs tab1;
    Out::os() << tab1 << "--------------------------------------------"
              << std::endl;
    Out::os() << tab1 << "k=" << k << std::endl;
    for (int i=0; i<sbasis.nterms(); i++)
    {
      Tabs tab2;
      Out::os() << tab2 << "\t";
      for (int j=0; j<sbasis.nterms(); j++)
      {
        Out::os() << sbasis.expectation(i,j,k) << " \t" ;
      }
      Out::os() << std::endl;
    }
  }


  Expr q = new SpectralExpr(sbasis, Q);
  Expr kappa = new SpectralExpr(sbasis, K);

  WatchFlag watch("watch");
  watch.setParam("integration setup", 0);
  watch.setParam("integration", 0);
  watch.setParam("fill", 0);
  watch.setParam("setup", 0);
  watch.setParam("evaluation", 0);
  watch.setParam("symbolic preprocessing", 0);
  watch.deactivate();

  /* Define the weak form */
  Expr eqn = Integral(interior, kappa*(dx*v)*(dx*u) + v*q, quad,watch);
  /* Define the Dirichlet BC */
  Expr bc = EssentialBC(leftPoint, v*u, quad);

  Out::os() << "done eq and bc " << std::endl;


  /* We can now set up the linear problem! */

  LinearProblem prob(mesh, eqn, bc, v, u, vecType); 


  LinearSolver<double> solver 
    = LinearSolverBuilder::createSolver("bicgstab.xml");

  Expr soln = prob.solve(solver);


  double totErrSq = 0.0;
  FieldWriter writer = new MatlabWriter("SpectralPoisson1DSoln");
  writer.addMesh(mesh);
  for (int i=0; i<sbasis.nterms(); i++)
  {
    Expr err = soln[i]-w[i];
    Expr errExpr = Integral(interior, err*err,
      new GaussianQuadrature(8));
    double errorSq = evaluateIntegral(mesh, errExpr);
    totErrSq += errorSq;
    Out::os() << "error norm [" << i << "] = " 
              << sqrt(errorSq) << std::endl << std::endl;
    writer.addField("u["+Teuchos::toString(i)+"]", 
      new ExprFieldWrapper(soln[i]));
  }
  writer.write();
      
        

  double tol = 1.0e-10;
  return SundanceGlobal::checkTest(sqrt(totErrSq), tol);
}

