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


#ifndef SUNDANCE_DIFFOPEVALUATOR_H
#define SUNDANCE_DIFFOPEVALUATOR_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Sundance 
{
class DiffOp;
class DiscreteFuncElementEvaluator;
    
/**
 *
 */
class DiffOpEvaluator : public UnaryEvaluator<DiffOp>
{
public:
  /** */
  DiffOpEvaluator(const DiffOp* expr,
    const EvalContext& context);

  /** */
  virtual ~DiffOpEvaluator(){;}

  /** */
  virtual void internalEval(const EvalManager& mgr,
    Array<double>& constantResults,
    Array<RCP<EvalVector> >& vectorResults) const ;

  /** We need a specialized resetting method for diff op
   * evaluators that also resets the discrete func evaluators
   * used in the functional chain rule */
  virtual void resetNumCalls() const ;

  /** */
  TEUCHOS_TIMER(evalTimer, "diff op evaluation");
private:

  Set<MultipleDeriv> increasedDerivs(const MultipleDeriv& mu,
    const Set<MultipleDeriv>& W1, int verb) const ;

  Set<MultipleDeriv> backedDerivs(const MultipleDeriv& mu,
    const Set<MultipleDeriv>& W1, int verb) const ;

  Deriv remainder(const MultipleDeriv& big, 
    const MultipleDeriv& little, int verb) const ;

  Array<int> isConstant_;

  Array<int> resultIndices_;
      
  Array<Array<int> > constantMonomials_;

  Array<Array<int> > vectorMonomials_;

  Array<Array<int> > constantFuncCoeffs_;

  Array<Array<int> > vectorFuncCoeffs_;

  Array<const DiscreteFuncElementEvaluator*> funcEvaluators_;

  /** Indices into the function evaluator table for the funcs
   * appearing with constant coeffs in the chain rule */
  Array<Array<int> > constantCoeffFuncIndices_;

  /** Indices into the list of multiindices for the funcs
   * appearing with constant coeffs in the chain rule */
  Array<Array<int> > constantCoeffFuncMi_;

  /** Indices into the function evaluator table for the funcs
   * appearing with vector coeffs in the chain rule */
  Array<Array<int> > vectorCoeffFuncIndices_;

  /** Indices into the list of multiindices for the funcs
   * appearing with vector coeffs in the chain rule */
  Array<Array<int> > vectorCoeffFuncMi_;
}; 
}


#endif
