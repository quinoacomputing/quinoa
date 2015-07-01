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

#ifndef SUNDANCE_SCALAREXPR_H
#define SUNDANCE_SCALAREXPR_H


#include "SundanceDefs.hpp"
#include "PlayaExceptions.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include "SundanceExprBase.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{

using namespace Sundance;
using namespace Teuchos;





/** */
class ScalarExpr : virtual public ExprBase
{
public:
  /** empty ctor */
  ScalarExpr();

  /** virtual destructor */
  virtual ~ScalarExpr() {;}


  /** Indicate whether this expression is constant in space */
  virtual bool isConstant() const {return false;}


  /** Indicate whether this expression is an immutable constant */
  virtual bool isImmutable() const {return false;}

  /** Indicate whether this expression is a "hungry"
   * differential operator that is awaiting an argument. */
  virtual bool isHungryDiffOp() const {return false;}

  /** Indicate whether the expression is independent of the given 
   * functions */
  virtual bool isIndependentOf(const Expr& u) const {return true;}

  /** 
   * Indicate whether the expression is nonlinear 
   * with respect to test functions */
  virtual bool isLinearInTests() const {return false;}

  /** 
   * Indicate whether every term in the expression contains test functions */
  virtual bool everyTermHasTestFunctions() const {return hasTestFunctions();}

  /** 
   * Indicate whether the expression contains test functions */
  virtual bool hasTestFunctions() const {return false;}
  /** 
   * Indicate whether the expression contains unknown functions */
  virtual bool hasUnkFunctions() const {return false;}

  /** Indicate whether the expression is linear in the given 
   * functions */
  virtual bool isLinearForm(const Expr& u) const {return false;}

  /** Indicate whether the expression is quadratic in the given 
   * functions */
  virtual bool isQuadraticForm(const Expr& u) const {return false;}

  /** Find all the unknown functions in this expression. */
  virtual void getUnknowns(Set<int>& unkID, Array<Expr>& unks) const {;}

  /** Find all the test functions in this expression. */
  virtual void getTests(Set<int>& varID, Array<Expr>& vars) const {;}

  /** Ordering operator for use in transforming exprs 
   * to standard form */
  virtual bool lessThan(const ScalarExpr* other) const = 0 ;


protected:
private:
};
}

#endif
