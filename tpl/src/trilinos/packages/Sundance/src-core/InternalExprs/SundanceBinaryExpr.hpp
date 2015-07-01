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

#ifndef SUNDANCE_BINARYEXPR_H
#define SUNDANCE_BINARYEXPR_H

#include "SundanceDefs.hpp"
#include "SundanceExprWithChildren.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * BinaryExpr is a base class for binary expressions, e.g., sums
 * and products. It provides a number of helper methods.
 */
class BinaryExpr : public ExprWithChildren
{
public:
  /** construct with left and right operands */
  BinaryExpr(const RCP<ScalarExpr>& left,
    const RCP<ScalarExpr>& right, int sign);

  /** virtual dtor */
  virtual ~BinaryExpr() {;}

  /** */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;

  /** */
  virtual XMLObject toXML() const ;

  /** */
  Expr left() const {return child(0);}

  /** */
  Expr right() const {return child(1);}

  /** */
  int sign() const {return sign_;}

  /** Downcast the left expr to an evaluatable expr */
  const EvaluatableExpr* leftEvaluatable() const 
    {return evaluatableChild(0);}

  /** Downcast the right expr to an evaluatable expr */
  const EvaluatableExpr* rightEvaluatable() const 
    {return evaluatableChild(1);}

  /** Downcast the left expr to a scalar expr */
  const ScalarExpr* leftScalar() const {return scalarChild(0);}

  /** Downcast the right expr to a scalar expr */
  const ScalarExpr* rightScalar() const {return scalarChild(1);}

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;

protected:

          

  /** */
  virtual bool parenthesizeSelf() const = 0 ;
  /** */
  virtual bool parenthesizeOperands() const = 0 ;
  /** */
  virtual const std::string& xmlTag() const = 0 ;
  /** */
  virtual const std::string& opChar() const = 0 ;



private:
          

  int sign_;
};
}

#endif
