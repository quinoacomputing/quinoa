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

#ifndef SUNDANCE_UNARYMINUS_H
#define SUNDANCE_UNARYMINUS_H

#include "SundanceDefs.hpp"
#include "SundanceUnaryExpr.hpp"
#include "SundanceUnaryMinusEvaluator.hpp"


namespace Sundance
{
  using namespace Sundance;
  using namespace Teuchos;
  

      /** */
      class UnaryMinus : public UnaryExpr,
                         GenericEvaluatorFactory<UnaryMinus, UnaryMinusEvaluator>
        {
        public:
          /** construct with the argument */
          UnaryMinus(const RCP<ScalarExpr>& arg);

          /** virtual dtor */
          virtual ~UnaryMinus() {;}


          /** */
          virtual std::ostream& toText(std::ostream& os, bool paren) const ;

          /** */
          virtual XMLObject toXML() const ;

          /** */
          virtual bool isLinear() const {return true;}

          /** 
           * Indicate whether the expression is linear 
           * with respect to test functions 
           */
          virtual bool isLinearInTests() const 
            {return evaluatableArg()->isLinearInTests();}
          

          /** Indicate whether the expression is linear in the given 
           * functions */
          virtual bool isLinearForm(const Expr& u) const 
            {
              return evaluatableArg()->isLinearForm(u);
            }

          /** Indicate whether the expression is at most 
           * quadratic in the given functions */
          virtual bool isQuadraticForm(const Expr& u) const
            {
              return evaluatableArg()->isQuadraticForm(u);
            }

          /** */
          virtual Set<MultiSet<int> > internalFindQ_W(int order, 
                                                   const EvalContext& context) const ;

          /** */
          virtual RCP<ExprBase> getRcp() {return rcp(this);}

        };
    }

#endif
