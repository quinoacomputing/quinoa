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

#ifndef SUNDANCE_CURVENORMEXPR_H
#define SUNDANCE_CURVENORMEXPR_H

#include "SundanceEvaluatorFactory.hpp"
#include "SundanceCurveNormEvaluator.hpp"

namespace Sundance
{
  using namespace Sundance;

  /** 
   * CurveNormExpr is an
   * expression that returns the Cartesian coordinate component
   * of the normal vector. Which dimensional component will be returned
   * can be specified in the constructor argument. <br>
   * The normal vector of the curve is always pointing outwards. This
   * outwards direction is in terms of the parameterized curve equation
   * from positive values to negative values.
   */
  class CurveNormExpr : public EvaluatableExpr,
                        public GenericEvaluatorFactory<CurveNormExpr, CurveNormEvaluator>
    {
    public:
      /** */
	  CurveNormExpr(int dir, const std::string& name="");

      /** */
      virtual ~CurveNormExpr() {;}

      /** */
      virtual XMLObject toXML() const ;

      /** */
      int dir() const {return dir_;}

      /** */
      const std::string& name() const {return name_;}

      /** Write a simple text description suitable 
       * for output to a terminal */
      virtual std::ostream& toText(std::ostream& os, bool paren) const 
        {os<<name(); return os;}
      
      /** */
      virtual Set<MultipleDeriv> 
      internalFindW(int order, const EvalContext& context) const ;
      
      /* we use the default implementation from the supper class */
      //virtual Set<MultipleDeriv>
      //internalFindC(int order, const EvalContext& context) const ;
      
      /*  we use the default implementation from the supper class */
      //virtual Set<MultipleDeriv>
      //internalFindV(int order, const EvalContext& context) const ;
      
      /** */
      virtual RCP<ExprBase> getRcp() {return rcp(this);}

      /** Ordering operator for use in transforming exprs to standard form */
      virtual bool lessThan(const ScalarExpr* other) const ;

      static std::string coordName(int dir, const std::string& name);

    private:
      int dir_;
      std::string name_;
    };
}

#endif
