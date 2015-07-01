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

#ifndef SUNDANCE_DISCRETEFUNCELEMENT_H
#define SUNDANCE_DISCRETEFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceFuncElementBase.hpp"
#include "SundanceDiscreteFuncEvaluator.hpp"
#include "SundanceDiscreteFuncDataStub.hpp"

namespace Sundance
{
using namespace Sundance;


using namespace Teuchos;




/** 
 * DiscreteFuncElement represents a scalar-valued element
 * of a (possibly) vector-valued discrete function. 
 *
 * DiscreteFuncElement is framework-independent. Any framework-specific
 * information should go in a subclass of DiscreteFuncDataStub.
 * The DiscreteFuncDataStub object can be accessed through the
 * <tt>master()</tt> method of this class.
 */
class DiscreteFuncElement : public virtual EvaluatableExpr,
                            public FuncElementBase,
                            public virtual GenericEvaluatorFactory<DiscreteFuncElement, DiscreteFuncElementEvaluator>
{
public:
  /** */
  DiscreteFuncElement(const RCP<DiscreteFuncDataStub>& data,
    const std::string& name,
    const std::string& suffix,
    const FunctionIdentifier& fid,
    int myIndexIntoVector);

  /** virtual destructor */
  virtual ~DiscreteFuncElement() {;}


  /** Get the data associated with the vector-valued function 
   * that contains this function element. */
  RCP<const DiscreteFuncDataStub> commonData() const {return commonData_;}

  /** Get the data associated with the vector-valued function 
   * that contains this function element. */
  DiscreteFuncDataStub* commonData() {return commonData_.get();}

  /** Get my index into the master's list of elements */
  int myIndex() const {return myIndex_;}

  /** Inform this function that it will need to be evaluated using the specified
   * multiIndex*/
  void addMultiIndex(const MultiIndex& newMi) const ;

  /**
   * Find the maximum differentiation order acting on discrete
   * functions in this expression. 
   */
  int maxDiffOrderOnDiscreteFunctions() const {return 0;}
      
  /**
   * Indicate whether this expression contains discrete functions.
   * This object is a discrete function, so return true.
   */
  virtual bool hasDiscreteFunctions() const {return true;}
      
  /**
   * Indicate whether this expression contains test functions.
   * This object is a discrete function, so return false.
   */
  virtual bool hasTestFunctions() const {return false;}

  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;
  /** */
  virtual Set<MultipleDeriv> 
  internalFindV(int order, const EvalContext& context) const ;
  /** */
  virtual Set<MultipleDeriv> 
  internalFindC(int order, const EvalContext& context) const ;

  /** */
  virtual RCP<Array<Set<MultipleDeriv> > > 
  internalDetermineR(const EvalContext& context,
    const Array<Set<MultipleDeriv> >& RInput) const ;

  /** */
  virtual XMLObject toXML() const ;

  /** */
  const Set<MultiIndex>& multiIndexSet() const {return miSet_;}

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** */
  bool lessThan(const ScalarExpr* other) const ;
      
private:

  RCP<DiscreteFuncDataStub> commonData_;

  mutable Set<MultiIndex> miSet_;

  int myIndex_;
      

};
}

#endif
