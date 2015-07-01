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

#ifndef SUNDANCE_UNKNOWNPARAMETERELEMENT_H
#define SUNDANCE_UNKNOWNPARAMETERELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceParameter.hpp"
#include "SundanceUnknownFuncDataStub.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"

namespace Sundance
{
using namespace Sundance;

using namespace Teuchos;




/** 
 * UnknownParameterElement represents an element of an unknown 
 * spatially-constant parameter
 */
class UnknownParameterElement : public UnknownFuncElement,
                                public SpatiallyConstantExpr
{
public:
  /** */
  UnknownParameterElement(const std::string& name,
    const std::string& suffix,
    const FunctionIdentifier& fid);


  /** virtual destructor */
  virtual ~UnknownParameterElement() {;}

      
  /** */
  Evaluator* createEvaluator(const EvaluatableExpr* expr,
    const EvalContext& context) const ;

      

  /** */
  void setValue(const double& value) {parameterValue()->setValue(value);}

  /** */
  const double& value() const {return parameterValue()->value();}

          
  /** */
  Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;
  /** */
  Set<MultipleDeriv> 
  internalFindC(int order, const EvalContext& context) const ;
  /** */
  Set<MultipleDeriv> 
  internalFindV(int order, const EvalContext& context) const ;

  /** */
  bool lessThan(const ScalarExpr* other) const ;

  /** */
  XMLObject toXML() const ;

  /** */
  bool isParameter() const {return true;}

  /** */
  RCP<ExprBase> getRcp() {return rcp(this);}
      
private:

  /** */
  const Parameter* parameterValue() const ;
  /** */
  Parameter* parameterValue() ;
      
};
}

#endif
