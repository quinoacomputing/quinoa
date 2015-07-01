
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

#ifndef SUNDANCE_SUMOFINTEGRALS_H
#define SUNDANCE_SUMOFINTEGRALS_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceOrderedHandle.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceMap.hpp"
#include "SundanceWatchFlag.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceParametrizedCurve.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;

using Sundance::Map;

class SpatiallyConstantExpr;

/** 
 * SumOfIntegrals represents a sum of integrals,
 * grouped by region and quadrature rule
 *
 * \f[
 * \sum_{d=0}^{N_d-1} \left[\sum_{q=0}^{N_{q,d}-1} 
 * \int_{\Omega_d,Q_{q,d}} g_{d,q}\right] 
 * \f] 
 *
 * 
 * 
 */
class SumOfIntegrals : public ScalarExpr
{
public:
  /** Construct given an integral over a single region */
  SumOfIntegrals(const RCP<CellFilterStub>& region,
    const Expr& expr,
    const RCP<QuadratureFamilyStub>& quad,
    const WatchFlag& watch);

  /** Construct given an integral over a single region */
  SumOfIntegrals(const RCP<CellFilterStub>& region,
    const Expr& expr,
    const RCP<QuadratureFamilyStub>& quad,
    const ParametrizedCurve& curve,
    const WatchFlag& watch);

  /** */
  virtual ~SumOfIntegrals(){;}

  /** Add another term to this integral */
  void addTerm(const RCP<CellFilterStub>& region,
    const Expr& expr,
    const RCP<QuadratureFamilyStub>& quad,
    const ParametrizedCurve& paramCurve,
    const WatchFlag& watch,
    int sign) ;

  /** Add this sum of integrals to another sum of integrals */
  void merge(const SumOfIntegrals* other, int sign) ;

  /** Multiply all terms in the sum by a constant */
  void multiplyByConstant(const SpatiallyConstantExpr* expr) ;

  /** Change the sign of all terms in the sum */
  void changeSign() ;

  /** Return the number of subregions */
  int numRQC() const {return rqcToExprMap_.size();}

  /** */
  const Sundance::Map<RegionQuadCombo, Expr>& rqcToExprMap() const
    {return rqcToExprMap_;}

  /** Return the set of unknown or variational
   * functions defined on region d */
  Set<int> funcsOnRegion(const OrderedHandle<CellFilterStub>& d, 
    const Set<int>& funcsSet) const ;

  /** Indicate whether the integral over 
   * region d contains any test functions */
  bool integralHasTestFunctions(const OrderedHandle<CellFilterStub>& d) const ;

  /** Return a null cell filter of a type consistent with the
   * other filters in this integral */
  RCP<CellFilterStub> nullRegion() const ;

  /** Indicate whether the expression is independent of the given 
   * functions */
  virtual bool isIndependentOf(const Expr& u) const ;

  /** Indicate whether the expression is linear in the given 
   * functions */
  virtual bool isLinearForm(const Expr& u) const ;

  /** Indicate whether the expression is quadratic in the given 
   * functions */
  virtual bool isQuadraticForm(const Expr& u) const ;

  /** 
   * Indicate whether the expression is nonlinear 
   * with respect to test functions */
  virtual bool isLinearInTests() const ;

  /** 
   * Indicate whether every term in the expression contains test functions */
  virtual bool everyTermHasTestFunctions() const ;

  /** 
   * Indicate whether the expression contains test functions */
  virtual bool hasTestFunctions() const ;
  

  /** Write a simple text description suitable 
   * for output to a terminal */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;


  /** Write in XML */
  virtual XMLObject toXML() const ;

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;

  /** Look for a term with an active watchpoint */
  bool hasWatchedTerm() const ;

  /** Find the maximum setup verbosity of the terms */
  int eqnSetSetupVerb() const ;
  

protected:
  /** */
  Expr filterSpectral(const Expr& ex) const ;
private:

  Sundance::Map<RegionQuadCombo, Expr> rqcToExprMap_;
};
}

#endif
