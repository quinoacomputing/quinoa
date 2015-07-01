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

#include "SundanceProblemTesting.hpp"

#include "SundanceMaximalCellFilter.hpp"
#include "SundanceDimensionalCellFilter.hpp"
#include "SundancePositionalCellPredicate.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshSource.hpp"
#include "SundanceMeshType.hpp"
#include "SundanceBasicSimplicialMeshType.hpp"
#include "SundancePartitionedLineMesher.hpp"
#include "SundancePartitionedRectangleMesher.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "SundanceGaussianQuadrature.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceLinearProblem.hpp"
#include "SundanceIntegral.hpp"

namespace Sundance
{

  using std::setprecision;
  using std::setw;

bool checkErrorNorms(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& numSoln,
  const Expr& exactSoln,
  const QuadratureFamily& quad,
  double L2Tol,
  double H1SemiTol,
  double H1Tol)
{
  Tabs tab;
  Tabs tab1;
  Expr df = numSoln - exactSoln;

  double L2Err = L2Norm(mesh, filter, df, quad);
  Out::root() << tab << "L2 norm check:" << std::endl
              << tab1 << setw(16) << setprecision(6) << L2Err 
              << " tol=" << setw(12) << L2Tol ;
  if (fabs(L2Err) > L2Tol) Out::root() << "   <==== FAIL!" << std::endl;
  else Out::root() << std::endl;

  double H1SemiErr = H1Seminorm(mesh, filter, df, quad);
  Out::root() << tab << "H1 seminorm check:" << std::endl
              << tab1 << setw(16) << setprecision(6) << H1SemiErr 
              << " tol=" << setw(12) << H1SemiTol ;
  if (fabs(H1SemiErr) > H1SemiTol) Out::root() << "   <==== FAIL!" << std::endl;
  else Out::root() << std::endl;


  double H1Err = H1Norm(mesh, filter, df, quad);
  Out::root() << tab << "H1 norm check:" << std::endl
              << tab1 << setw(16) << setprecision(6) << H1Err 
              << " tol=" << setw(12) << H1Tol;
  if (fabs(H1Err) > H1Tol) Out::root() << "   <==== FAIL!" << std::endl;
  else Out::root() << std::endl;

  return (fabs(L2Err) <= L2Tol) 
    && (fabs(H1SemiErr) <= H1SemiTol) 
    && (fabs(H1Err) <= H1Tol) ;
}

double fitPower(const Array<double>& h, const Array<double>& err)
{
  Array<double> x(h.size());
  Array<double> y(h.size());
  double xBar = 0.0;
  double yBar = 0.0;
  for (int i=0; i<h.size(); i++)
  {
    x[i] = log(h[i]);
    y[i] = log(err[i]);
    xBar += x[i];
    yBar += y[i];
  }

  xBar /= h.size();
  yBar /= h.size();
  
  double u = 0.0;
  double v = 0.0;
  for (int i=0; i<h.size(); i++)
  {
    u += (x[i]-xBar)*(y[i]-yBar);
    v += pow(x[i]-xBar,2.0);
  }

  return u/v;
}

/* -------- LineDomain ----------------- */

LineDomain::LineDomain(const Array<int>& nx)
  : a_(0.0), b_(1.0), nx_(nx), interior_(new MaximalCellFilter()),
    left_(), right_(), mesh_()
{init();}


LineDomain::LineDomain(double a, double b, const Array<int>& nx)
  : a_(a), b_(b), nx_(nx), interior_(new MaximalCellFilter()),
    left_(), right_(), mesh_()
{init();}

void LineDomain::init()
{
  int np = MPIComm::world().getNProc();
  MeshType meshType = new BasicSimplicialMeshType();
  mesh_.resize(nx_.size());
  for (int i=0; i<nx_.size(); i++)
  {
    MeshSource mesher = new PartitionedLineMesher(a_, b_, np*nx_[i], meshType);
    mesh_[i] = mesher.getMesh();
  }
  
  CellFilter points = new DimensionalCellFilter(0);
  left_ = points.subset(new CoordinateValueCellPredicate(0,a_));
  right_ = points.subset(new CoordinateValueCellPredicate(0,b_));
}

/* -------- RectangleDomain ----------------- */

RectangleDomain::RectangleDomain(const Array<int>& n)
  : ax_(0.0), bx_(1.0), nx_(n), 
    ay_(0.0), by_(1.0), ny_(n), 
    interior_(new MaximalCellFilter()),
    north_(), south_(), east_(), west_(),
    mesh_()
{init();}

RectangleDomain::RectangleDomain(
  double ax, double bx, const Array<int>& nx,
  double ay, double by, const Array<int>& ny
  )
  : ax_(ax), bx_(bx), nx_(nx), 
    ay_(ay), by_(by), ny_(ny), 
    interior_(new MaximalCellFilter()),
    north_(), south_(), east_(), west_(),
    mesh_()
{init();}


void RectangleDomain::init()
{
  TEUCHOS_TEST_FOR_EXCEPT(nx_.size() != ny_.size());
  int np = MPIComm::world().getNProc();
  int npx = -1;
  int npy = -1;

  PartitionedRectangleMesher::balanceXY(np, &npx, &npy);
  TEUCHOS_TEST_FOR_EXCEPT(npx < 1);
  TEUCHOS_TEST_FOR_EXCEPT(npy < 1);
  TEUCHOS_TEST_FOR_EXCEPT(npx * npy != np);

  MeshType meshType = new BasicSimplicialMeshType();
  mesh_.resize(nx_.size());

  for (int i=0; i<nx_.size(); i++)
  {
    MeshSource mesher 
      = new PartitionedRectangleMesher(
        ax_, bx_, nx_[i], npx,
        ay_, by_, ny_[i], npy,
        meshType);
    mesh_[i] = mesher.getMesh();
  }
  
  CellFilter edges = new DimensionalCellFilter(1);
  north_ = edges.subset(new CoordinateValueCellPredicate(1,by_));
  south_ = edges.subset(new CoordinateValueCellPredicate(1,ay_));
  west_ = edges.subset(new CoordinateValueCellPredicate(0,ax_));
  east_ = edges.subset(new CoordinateValueCellPredicate(0,bx_));
}





/* -------- norm calculators ----------------- */

Array<double> L2NormCalculator::computeNorms(
  const ForwardProblemTestBase* prob,
  int meshIndex,
  const Expr& numSoln, const Expr& exactSoln) const
{
  Expr errFunc = (numSoln - exactSoln).flatten();

  Mesh mesh = prob->getMesh(meshIndex);
  CellFilter interior = prob->interior();

  Array<int> p = prob->pExpected();
  TEUCHOS_TEST_FOR_EXCEPTION(p.size() != errFunc.size(),
    std::runtime_error,
    "size mismatch between array of expected orders (p=" << p << ") and "
    "array of solutions: " << errFunc.size());
 
  Array<double> rtn(p.size());
  for (int i=0; i<errFunc.size(); i++)
  {
    QuadratureFamily quad = new GaussianQuadrature(2*p[i]);
    double L2Err = L2Norm(mesh, interior, errFunc[i], quad);
    rtn[i] = L2Err;
  }  

  return rtn;
}
  

Array<LPTestSpec> LPTestBase::specs() const
{
  return tuple(
    LPTestSpec("amesos.xml", 1.0e-10, makeSet<int>(1)),
    LPTestSpec("aztec-ifpack.xml", 1.0e-10),
    LPTestSpec("aztec-ml.xml", 1.0e-10),
    LPTestSpec("belos-ifpack.xml", 1.0e-8),
    LPTestSpec("belos-ml.xml", 1.0e-10)
    );
}

/** \relates LPTestSpec */
std::ostream& operator<<(std::ostream& os, const LPTestSpec& spec)
{
  os << "LPTestSpec(tol=" << spec.tol() << ", solver=" << spec.solverFile()
     << std::endl;
  return os;
}




bool LPTestSuite::run() const 
{
  int np = MPIComm::world().getNProc();

  bool allOK = true;

  for (int i=0; i<tests_.size(); i++)
  {
    Tabs tab(0);
    Array<LPTestSpec> specs = tests_[i]->specs();
    for (int j=0; j<specs.size(); j++)
    {
      if (specs[j].nProcIsAllowed(np))
      {
        Out::root() << std::endl;
        Out::root() << std::endl;
        Out::root() << std::endl;
        Out::root() << tab 
                    << "-------------------------------------" 
          "-------------------------------------" 
                    << std::endl;
        Out::root() << tab << "running test " << tests_[i]->name()
                    << " with spec " << specs[j] << std::endl;
        Out::root() << tab 
                    << "-------------------------------------" 
          "-------------------------------------" 
                    << std::endl;

        std::string solverFile = specs[j].solverFile();
        double tol = specs[j].tol();
        bool pass = tests_[i]->run(solverFile, tol);
        allOK = pass && allOK;
      }
      else
      {
        Out::root() << tab << "skipping test " << tests_[i]->name()
                    << " with spec=" << specs[j] << std::endl;
      }
    }
    Out::root() << std::endl;
    Out::root() << std::endl;
  }
  return allOK;
}


LPTestSuite::LPTestSuite()
  : tests_(), testSpecs_() {}


void LPTestSuite::registerTest(const RCP<LPTestBase>& test) 
{
  tests_.append(test);
}

VectorType<double> ForwardProblemTestBase::vecType() const
{
  return new EpetraVectorType();
}

Expr ForwardProblemTestBase::coord(int d) const 
{
  TEUCHOS_TEST_FOR_EXCEPT(d<0 || d>2);
  return new CoordExpr(d);
}

double ForwardProblemTestBase::cellSize(int i) const 
{
  Expr h = new CellDiameterExpr();
  Expr hExpr = Integral(interior(), h, new GaussianQuadrature(1));
  Expr AExpr = Integral(interior(), 1.0, new GaussianQuadrature(1));
  
  double area = evaluateIntegral(getMesh(i), AExpr);
  double hMean = evaluateIntegral(getMesh(i), hExpr)/area;

  return hMean;
}

RCP<ErrNormCalculatorBase> ForwardProblemTestBase::normCalculator() const
{
  return rcp(new L2NormCalculator());

}
bool ForwardProblemTestBase::run(const std::string& solverFile, 
  double tol) const
{
  if (numMeshes()==1) return runSingleTest(solverFile, tol);
  else return runTestSequence(solverFile, tol);
}


bool ForwardProblemTestBase::runTestSequence(const std::string& solverFile, 
  const double& pTol) const
{
  Tabs tab(0);
  Out::root() << tab << "Running test sequence " << name() << std::endl;

  LinearSolver<double> solver 
    = LinearSolverBuilder::createSolver(solverFile);

  Expr exact = exactSoln();

  Array<double> hList(numMeshes());
  Array<Array<double> > errList(numMeshes());

  for (int i=0; i<numMeshes(); i++)
  {
    Tabs tab1;
    Out::root() << tab1 << "running on mesh #" << i << std::endl;
    Expr soln;
    bool solveOK = solve(getMesh(i), solver, soln);
    if (!solveOK) return false;

    postRunCallback(i, getMesh(i), solverFile, soln);

    double h = cellSize(i);
    Array<double> err = normCalculator()->computeNorms(
      this, i, soln, exact);
    hList[i] = h;
    errList[i] = err;
  }

  bool allPass = true;
  Out::root() << tab << "Error results: " << std::endl;
  for (int f=0; f < pExpected().size(); f++)
  {
    Tabs tab1;
    Out::root() << tab1 << "Variable #" << f << std::endl; 
    Out::root() << std::endl << tab1 << "Error norm versus h" << std::endl;
    Array<double> err;
    for (int i=0; i<numMeshes(); i++)
    {
      Tabs tab2;
      Out::root() << tab2 << setw(20) << setprecision(5) << hList[i] 
                  << " " << setw(20) << setprecision(5) << errList[i][f] 
                  << std::endl;
      err.append(errList[i][f]); 
    }
    
    double p = fitPower(hList, err);
    Out::root() << tab << "Measured exponent: " << p << std::endl;
    Out::root() << tab << "Expected exponent: " << pExpected()[f] << std::endl;
    double pErr = fabs(p - pExpected()[f]) ;
    Out::root() << tab << "Difference: " << setw(12) << pErr
                << " Tolerance: " << setw(12) << pTol;
    bool pass = pErr <= pTol;
    if (!pass) Out::root() << "  <==== FAIL!";
    Out::root() << std::endl;
    allPass = allPass && pass;
  }

  return allPass;
}



bool ForwardProblemTestBase::runSingleTest(const std::string& solverFile, 
  const double& tol) const
{
  Tabs tab(0);

  LinearSolver<double> solver 
    = LinearSolverBuilder::createSolver(solverFile);

  Expr exact = exactSoln();

  Expr soln;
  bool solveOK = solve(getMesh(0), solver, soln);
  if (!solveOK) return false;
  
  Array<double> err = normCalculator()->computeNorms(this, 0, 
    soln, exact);
  
  for (int i=0; i<err.size(); i++)
  {
    Out::root() << tab << setw(20) << setprecision(5) << err[i] 
                << " " << setw(20) << setprecision(5) << tol 
              << std::endl;
  }
  return err[0] <= tol;
}

bool LPTestBase::solve(const Mesh& mesh,
  const LinearSolver<double>& solver,
  Expr& soln) const
{
  Tabs tab(0);
  LinearProblem::solveFailureIsFatal() = false;
  SolverState<double> state = prob(mesh).solve(solver, soln);

  bool converged = (state.finalState() == SolveConverged) ;
  if (!converged)
  {
    Out::root() << tab << "solve failed to converge!" << std::endl;
    Tabs tab1;
    Out::root() << tab1 << "state = " << state << std::endl;
  }
  return converged;
}

LP1DTestBase::LP1DTestBase(const Array<int>& nx)
  : domain_(nx) {}

LP1DTestBase::LP1DTestBase(double a, double b, const Array<int>& nx)
  : domain_(a, b, nx) {}

LPRectTestBase::LPRectTestBase(const Array<int>& n)
  : domain_(n) {}
}
