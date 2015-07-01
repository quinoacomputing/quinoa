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

#ifndef SUNDANCE_FUNCTIONALPOLYNOMIAL_H
#define SUNDANCE_FUNCTIONALPOLYNOMIAL_H

#include "SundanceDefs.hpp"
#include "SundanceEvaluatableExpr.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




/**
 * Specialized class for representing polynomials in symbolic
 * functions and their derivatives. 
 */
class FunctionalPolynomial : public EvaluatableExpr
{
public:
  /** ctor */
  FunctionalPolynomial(const RCP<ScalarExpr>& expr);
  /** ctor */
  FunctionalPolynomial(const Map<int, RCP<ScalarExpr> >& funcs,
    const Map<int, Set<MultiIndex> >& funcMultiIndices,
    const Array<Map<MultipleDeriv, RCP<ScalarExpr> > > & coeffs);

  /** virtual destructor */
  virtual ~FunctionalPolynomial() {;}


  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;

      

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** */
  virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;

  /** */
  RCP<FunctionalPolynomial> addPoly(const FunctionalPolynomial* other,
    int sign) const ;

  /** */
  RCP<FunctionalPolynomial> multiplyPoly(const FunctionalPolynomial* other) const ;

  /** */
  RCP<FunctionalPolynomial> multiplyScalar(const RCP<ScalarExpr>& alpha) const ;

  /** */
  RCP<FunctionalPolynomial> addFunction(const RCP<ScalarExpr>& u,
    int sign) const ;

  /** */
  static bool isConvertibleToPoly(const ScalarExpr* expr) ;

  /** */
  static RCP<FunctionalPolynomial> toPoly(const RCP<ScalarExpr>& expr);


  /** Write a simple text description suitable 
   * for output to a terminal */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;
      
  /** Write in XML */
  virtual XMLObject toXML() const ;

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;

private:

  /** */
  Map<int, RCP<ScalarExpr> > funcs_;

  /** */
  Map<int, Set<MultiIndex> > funcMultiIndices_;

  /** */
  Array<Map<MultipleDeriv, RCP<ScalarExpr> > > coeffs_;

  /** */
  Array<Set<MultipleDeriv> > keys_;

  /** */
  Set<Deriv> findFuncsForSummation(const Set<MultipleDeriv>& prevSet,
    const MultipleDeriv& thisSet) const ;
      
  /**    
   * Given a term's key, find the term that will appear as
   * its successor in the evaluation recurrence.
   */
  MultipleDeriv successorTerm(const MultipleDeriv& md) const ;

  /** */
  void stepRecurrence(int level, const Map<MultipleDeriv, std::string>& sPrev,
    Map<MultipleDeriv, std::string>& sCurr) const ;

  /** */
  std::string evalString() const ;
};
}

#endif
