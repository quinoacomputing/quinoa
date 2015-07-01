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

#ifndef SUNDANCE_USERDEFOPELEMENT_H
#define SUNDANCE_USERDEFOPELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceMultiSet.hpp"
#include "SundanceSet.hpp"
#include "SundanceMultipleDeriv.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceUserDefFunctorElement.hpp"
#include "SundanceUserDefOpEvaluator.hpp"


namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;





class UserDefOp;

/** 
 * Scalar element of a vector-valued user-defined expression.
 */
class UserDefOpElement : virtual public ExprWithChildren
{
public:
  /** */
  UserDefOpElement(const Array<RCP<ScalarExpr> >& args,
    const RCP<Sundance::Map<EvalContext, RCP<const UserDefOpCommonEvaluator> > >& evalMap,
    const RCP<const UserDefFunctorElement>& functorElement);

  /** virtual destructor */
  virtual ~UserDefOpElement() {;}

  /** Return the index of this element into 
   * the list-valued user defined op */
  int myIndex() const {return functorElement_->myIndex();}


  /** Write self in text form */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;

  /** Write in XML */
  virtual XMLObject toXML() const ;

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** Access to the functor underlying this object */
  const UserDefFunctorElement* functorElement() const 
    {return functorElement_.get();}

  /** */
  void reset() const ;

  /** */
  Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;

  /** */
  virtual void getArgDerivIndices(const Array<int>& orders,
    Sundance::Map<MultiSet<int>, int>& varArgDerivs,
    Sundance::Map<MultiSet<int>, int>& constArgDerivs) const ;

protected:
  /** Get an evaluator that will be common to all vector elements of 
   * this operator */
  RCP<const UserDefOpCommonEvaluator> 
  getCommonEvaluator(const EvalContext& context) const ;
          
private:
  mutable RCP<Sundance::Map<EvalContext, RCP<const UserDefOpCommonEvaluator> > > commonEvaluatorsMap_;
  const RCP<const UserDefFunctorElement> functorElement_;
};
}

#endif
