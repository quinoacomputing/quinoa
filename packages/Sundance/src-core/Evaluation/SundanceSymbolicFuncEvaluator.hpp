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

#ifndef SUNDANCE_SYMBOLICFUNCEVALUATOR_H
#define SUNDANCE_SYMBOLICFUNCEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceSubtypeEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"


namespace Sundance 
{
class Parameter;

class DiscreteFuncElementEvaluator; 
class SymbolicFuncElement;
class ConstantEvaluator; 

/** 
 *
 */
class SymbolicFuncElementEvaluator 
  : public SubtypeEvaluator<SymbolicFuncElement>
{
public:
  /** */
  SymbolicFuncElementEvaluator(const SymbolicFuncElement* expr, 
    const EvalContext& context);

  /** */
  virtual ~SymbolicFuncElementEvaluator(){;}

  /** */
  virtual void internalEval(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const ;

  /** */
  TEUCHOS_TIMER(symbolicFuncEvalTimer, "symbolic function evaluation");

  /** */
  const DiscreteFuncElementEvaluator* dfEval() const {return dfEval_;}
  /** */
  const ConstantEvaluator* pEval() const {return pEval_;}


  /** Reset the number of calls to zero. This should be called
   * at the beginning of every new evaluation cycle. */
  virtual void resetNumCalls() const ;
private:
  Array<MultiIndex> mi_;
  Array<int> spatialDerivPtrs_;
  Array<int> onePtrs_;
  Array<int> paramValuePtrs_;
  const DiscreteFuncElement* df_;
  const Parameter* p_;
  const DiscreteFuncElementEvaluator* dfEval_;
  const ConstantEvaluator* pEval_;
  Array<string> stringReps_;
};

    

}

#endif
