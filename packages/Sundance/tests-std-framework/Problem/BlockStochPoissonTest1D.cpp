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
#include "SundanceStochBlockJacobiSolver.hpp"
#ifdef HAVE_SUNDANCE_STOKHOS
#include "Stokhos_HermiteBasis.hpp"
#endif



bool BlockStochPoissonTest1D()
{
  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Read a mesh */
  MeshType meshType = new BasicSimplicialMeshType();
  int nx = 32;
  MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, 
    meshType);
  Mesh mesh = mesher.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
  CellFilter pts = new DimensionalCellFilter(0);
  CellFilter left = pts.subset(new CoordinateValueCellPredicate(0,0.0));
  CellFilter right = pts.subset(new CoordinateValueCellPredicate(0,1.0));

  Expr x = new CoordExpr(0);

  /* Create the stochastic coefficients */
  int nDim = 1;
  int order = 6;
#ifdef HAVE_SUNDANCE_STOKHOS
  Out::root() << "using Stokhos hermite basis" << std::endl;
  SpectralBasis pcBasis = new Stokhos::HermiteBasis<int,double>(order);
#else
  Out::root() << "using George's hermite basis" << std::endl;
  SpectralBasis pcBasis = new HermiteSpectralBasis(nDim, order);
#endif
    
  Array<Expr> q(pcBasis.nterms());
  Array<Expr> kappa(pcBasis.nterms());
  Array<Expr> uEx(pcBasis.nterms());

  double a = 0.1;

  q[0] = -2 + pow(a,2)*(4 - 9*x)*x - 2*pow(a,3)*(-1 + x)*(1 + 3*x*(-3 + 4*x));
  q[1] = -(a*(-3 + 10*x + 2*a*(-1 + x*(8 - 9*x +
          a*(-4 + 3*(5 - 4*x)*x + 12*a*(-1 + x)*(1 + 5*(-1 + x)*x))))));
  q[2] = a*(-4 + 6*x + a*(1 - x*(2 + 3*x) + a*(4 - 28*x + 30*pow(x,2))));
  q[3] = -(pow(a,2)*(-3 + x*(20 - 21*x +
        a*(-4 + 3*(5 - 4*x)*x + 24*a*(-1 + x)*(1 + 5*(-1 + x)*x)))));
  q[4] = pow(a,3)*(1 + x*(-6 + x*(3 + 4*x)));
  q[5] = -4*pow(a,4)*(-1 + x)*x*(1 + 5*(-1 + x)*x);
  q[6] = 0.0;

  uEx[0] = -((-1 + x)*x);
  uEx[1] = -(a*(-1 + x)*pow(x,2));
  uEx[2] = a*pow(-1 + x,2)*x;
  uEx[3] = pow(a,2)*pow(-1 + x,2)*pow(x,2);
  uEx[4] = 0.0;
  uEx[5] = 0.0;
  uEx[6] = 0.0;

  kappa[0] = 1.0;
  kappa[1] = a*x;
  kappa[2] = -(pow(a,2)*(-1 + x)*x);

  kappa[3] = 1.0; // unused
  kappa[4] = 1.0; // unused
  kappa[5] = 1.0; // unused
  kappa[6] = 1.0; // unused


  Array<Expr> uBC(pcBasis.nterms());
  for (int i=0; i<pcBasis.nterms(); i++) uBC[i] = 0.0;

  int L = nDim+2;
  int P = pcBasis.nterms();
  Out::os() << "L = " << L << std::endl;
  Out::os() << "P = " << P << std::endl;
    
  /* Create the unknown and test functions. Do NOT use the spectral
   * basis here */
  Expr u = new UnknownFunction(new Lagrange(4), "u");
  Expr v = new TestFunction(new Lagrange(4), "v");

  /* Create differential operator and coordinate function */
  Expr dx = new Derivative(0);
  Expr grad = dx;


  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad = new GaussianQuadrature(12);

  /* Now we create problem objects to build each $K_j$ and $f_j$.
   * There will be L matrix-vector pairs */
  Array<Expr> eqn(P);
  Array<Expr> bc(P);
  Array<LinearProblem> prob(P);
  Array<LinearOperator<double> > KBlock(L);
  Array<Vector<double> > fBlock(P);
  Array<Vector<double> > solnBlock;

  for (int j=0; j<P; j++)
  {
    eqn[j] = Integral(interior, kappa[j]*(grad*v)*(grad*u) + v*q[j], quad);
    bc[j] = EssentialBC(left+right, v*(u-uBC[j]), quad);
    prob[j] = LinearProblem(mesh, eqn[j], bc[j], v, u, vecType); 
    if (j<L) KBlock[j] = prob[j].getOperator();
    fBlock[j] = -1.0*prob[j].getSingleRHS();
  }

  /* Read the solver to be used on the diagonal blocks */
  LinearSolver<double> diagSolver 
    = LinearSolverBuilder::createSolver("amesos.xml");

    
  double convTol = 1.0e-12;
  int maxIters = 30;
  int verb = 1;
  StochBlockJacobiSolver solver(diagSolver, pcBasis,
    convTol, maxIters, verb);
    
  solver.solve(KBlock, fBlock, solnBlock);

  /* write the solution */
  FieldWriter w = new MatlabWriter("Stoch1D");
  w.addMesh(mesh);
  DiscreteSpace discSpace(mesh, new Lagrange(4), vecType);
  for (int i=0; i<P; i++)
  {
    L2Projector proj(discSpace, uEx[i]);
    Expr ue_i = proj.project();
    Expr df = new DiscreteFunction(discSpace, solnBlock[i]);
    w.addField("u["+ Teuchos::toString(i)+"]", 
      new ExprFieldWrapper(df));
    w.addField("uEx["+ Teuchos::toString(i)+"]", 
      new ExprFieldWrapper(ue_i));
  }
  w.write();

  double totalErr2 = 0.0;
  DiscreteSpace discSpace4(mesh, new Lagrange(4), vecType);
  for (int i=0; i<P; i++)
  {
    Expr df = new DiscreteFunction(discSpace4, solnBlock[i]);
    Expr errExpr = Integral(interior, pow(uEx[i]-df, 2.0), quad);
    Expr scaleExpr = Integral(interior, pow(uEx[i], 2.0), quad);
    double errSq = evaluateIntegral(mesh, errExpr);
    double scale = evaluateIntegral(mesh, scaleExpr);
    if (scale > 0.0) 
      Out::os() << "mode i=" << i << " error=" << sqrt(errSq/scale) << std::endl;
    else
      Out::os() << "mode i=" << i << " error=" << sqrt(errSq) << std::endl;
  }
    
  double tol = 1.0e-12;
    
  return SundanceGlobal::checkTest(sqrt(totalErr2), tol);
}

    

    

    
    

    
    
