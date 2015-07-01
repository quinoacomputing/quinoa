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

#ifndef SUNDANCE_PROBLEMTESTING_H
#define SUNDANCE_PROBLEMTESTING_H

#include "SundanceFunctional.hpp"
#include "SundanceLinearProblem.hpp"

namespace Sundance
{

using namespace Teuchos;


/** 
 * This function checks the L2 and H1 norms and the H1 seminorm of
 * an error against a specified tolerance.
 */
bool checkErrorNorms(
  const Mesh& mesh,
  const CellFilter& filter,
  const Expr& numSoln,
  const Expr& exactSoln,
  const QuadratureFamily& quad,
  double L2Tol,
  double H1SemiTol,
  double H1Tol);

/** 
 * Compute the exponent \f$p\f$ that best fits the measured errors to a 
 * power law \f$\epsilon=A h^p\f$. 
 */
double fitPower(const Array<double>& h, const Array<double>& err);


/** 
 * This class bundles together a sequence of uniform meshes
 * of the 1D interval [a,b] with cell
 * filters defining the interior and boundaries. It is intended for quick 
 * and reliable setup of 1D test problems.
 */
class LineDomain
{
public:
  /** */
  LineDomain(const Array<int>& nx);

  /** */
  LineDomain(double a, double b, const Array<int>& nx);

  /** */
  int numMeshes() const {return mesh_.size();}

  /** */
  const CellFilter& left() const {return left_;}

  /** */
  const CellFilter& right() const {return right_;}

  /** */
  const CellFilter& interior() const {return interior_;}

  /** */
  const Mesh& mesh(int i) const {return mesh_[i];}

  /** */
  double a() const {return a_;}

  /** */
  double b() const {return b_;}
  
  /** */
  int nx(int i) const {return nx_[i];}

private:
  void init();

  double a_;
  double b_;
  Array<int> nx_;
  CellFilter interior_;
  CellFilter left_;
  CellFilter right_;
  Array<Mesh> mesh_;
};


/** 
 * This class bundles together a sequence of uniform meshes
 * of the 2D rectangle [ax,bx] times [ay,by] with cell
 * filters defining the interior and boundaries. It is intended for quick 
 * and reliable setup of 2D test problems.
 */
class RectangleDomain
{
public:
  /** Contruct the unit square */
  RectangleDomain(const Array<int>& n);

  /** Construct */
  RectangleDomain(
    double ax, double bx, const Array<int>& nx,
    double ay, double by, const Array<int>& ny
    );

  /** */
  int numMeshes() const {return mesh_.size();}

  /** */
  const CellFilter& north() const {return north_;}

  /** */
  const CellFilter& south() const {return south_;}

  /** */
  const CellFilter& east() const {return east_;}

  /** */
  const CellFilter& west() const {return west_;}

  /** */
  const CellFilter& interior() const {return interior_;}

  /** */
  const Mesh& mesh(int i) const {return mesh_[i];}

  /** */
  double ax() const {return ax_;}

  /** */
  double bx() const {return bx_;}
  
  /** */
  int nx(int i) const {return nx_[i];}

  /** */
  double ay() const {return ay_;}

  /** */
  double by() const {return by_;}
  
  /** */
  int ny(int i) const {return ny_[i];}

private:
  void init();

  double ax_;
  double bx_;
  Array<int> nx_;
  double ay_;
  double by_;
  Array<int> ny_;
  CellFilter interior_;
  CellFilter north_;
  CellFilter south_;
  CellFilter east_;
  CellFilter west_;
  Array<Mesh> mesh_;
};



/** */
class LPTestSpec
{
public:
  /** */
  LPTestSpec() {;}
  /** */
  LPTestSpec(const std::string& solverFile, double tol)
    : hasProcRestriction_(false), allowedProcNumbers_(),
      solverFile_(solverFile), tol_(tol){}

  /** */
  LPTestSpec(const std::string& solverFile, double tol, 
    const Set<int>& allowedProcs) 
    : hasProcRestriction_(true), allowedProcNumbers_(allowedProcs),
      solverFile_(solverFile), tol_(tol){}

  /** */
  const double& tol() const {return tol_;}

  /** */
  const std::string& solverFile() const {return solverFile_;}

  /** */
  bool nProcIsAllowed(int np) const
    {
      if (!hasProcRestriction_) return true;
      return allowedProcNumbers_.contains(np);
    }

private:
  bool hasProcRestriction_;

  Set<int> allowedProcNumbers_;

  std::string solverFile_;

  double tol_;
};


/** \relates LPTestSpec */
std::ostream& operator<<(std::ostream& os, const LPTestSpec& spec);

class ForwardProblemTestBase;

/**
 * Base class for objects that compute error norms. 
 */
class ErrNormCalculatorBase
{
public:
  /** Compute the error norm. In vector-valued problems we may need to
   * compute multiple norms, so the return type is an array. */
  virtual Array<double> computeNorms(const ForwardProblemTestBase* prob,
    int meshIndex,
    const Expr& numSoln, const Expr& exactSoln) const = 0 ;
};

/**
 * Object to compute L2 norms of errors 
 */
class L2NormCalculator : public ErrNormCalculatorBase
{
public:
  /** */
  L2NormCalculator() {}

  /** */
  virtual Array<double> computeNorms(const ForwardProblemTestBase* prob,
    int meshIndex,
    const Expr& numSoln, const Expr& exactSoln) const ;

};

/** 
 * 
 */
class ForwardProblemTestBase
{
public:
  /** */
  virtual bool run(const std::string& solverFile, double tol) const ;

  /** */
  virtual std::string name() const = 0 ;

  /** */
  virtual Expr exactSoln() const = 0 ;

  /** */
  virtual VectorType<double> vecType() const ;

  /** */
  virtual Expr coord(int d) const ;

  /** */
  virtual Mesh getMesh(int i) const = 0 ;

  /** */
  virtual CellFilter interior() const = 0 ;

  /** */
  virtual RCP<ErrNormCalculatorBase> normCalculator() const ;

  /** 
   * Solve the problem on the \f$i\f$-th mesh. Return a bool indicating whether
   * the solve succeeded. 
   */
  virtual bool solve(const Mesh& mesh, const LinearSolver<double>& solver,
    Expr& soln) const = 0 ;

  /** */
  virtual int numMeshes() const = 0 ;

  /** 
   * Return the average cell size on the \f$i\f$-th mesh.
   */
  virtual double cellSize(int i) const ;

  /** 
   * Return the order of accuracy expected for the solution. If the problem
   * is vector-valued, an array of expected orders is returned.  
   */
  virtual Array<int> pExpected() const = 0 ;

  /** 
   * Hook for a user-defined operation after each run. This can be used,
   * for example, to write viz output. Default is a no-op.
   */
  virtual void postRunCallback(int meshID, const Mesh& mesh,
    const string& solverFile,
    const Expr& soln) const {} 

private:
  /** */
  bool runSingleTest(const std::string& solverFile, const double& tol) const ;

  /** */
  bool runTestSequence(const std::string& solverFile, const double& tol) const ;
};

/** */
class LPTestBase : public ForwardProblemTestBase
{
public:

  /** */
  virtual Array<LPTestSpec> specs() const ;

  /** */
  virtual LinearProblem prob(const Mesh& mesh) const = 0 ;

  /** */
  virtual bool solve(const Mesh& mesh, 
    const LinearSolver<double>& solver,
    Expr& soln) const ;
};



/** */
class LP1DTestBase : public LPTestBase
{
public:
  /** */
  LP1DTestBase(const Array<int>& nx);

  /** */
  LP1DTestBase(double a, double b, const Array<int>& nx);

  /** */
  CellFilter interior() const {return domain_.interior();}

  /** */
  Mesh getMesh(int i) const {return domain_.mesh(i);}

  /** */
  const LineDomain& domain() const {return domain_;}

  /** */
  int numMeshes() const {return domain_.numMeshes();}
  
private:
  LineDomain domain_;
};


/** */
class LPRectTestBase : public LPTestBase
{
public:
  /** */
  LPRectTestBase(const Array<int>& n);

  /** */
  CellFilter interior() const {return domain_.interior();}

  /** */
  Mesh getMesh(int i) const {return domain_.mesh(i);}

  /** */
  const RectangleDomain& domain() const {return domain_;}

  /** */
  int numMeshes() const {return domain_.numMeshes();}
  
private:
  RectangleDomain domain_;
};


/** */
class LPTestSuite
{
public:
  /** */
  LPTestSuite();

  /** */
  void registerTest(const RCP<LPTestBase>& test) ;

  /** */
  bool run() const ;
  
private:
  Array<RCP<LPTestBase> > tests_;
  Array<Array<LPTestSpec> > testSpecs_;
};


}

#endif

