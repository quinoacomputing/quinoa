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

#include "SundanceDerivOfSymbFunc.hpp"

#include "SundanceTestFuncElement.hpp"

#include "SundanceCoordExpr.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


DerivOfSymbFunc::DerivOfSymbFunc(const MultiIndex& op, 
  const RCP<ScalarExpr>& arg)
  : DiffOp(op, arg), argFid_(FunctionIdentifier())
{
  const SymbolicFuncElement* f 
    = dynamic_cast<const SymbolicFuncElement*>(evaluatableArg());
  TEUCHOS_TEST_FOR_EXCEPTION(f==0, std::logic_error, "argument to DerivOfSymbFunc ctor "
                     "is not a symbolic function");
  argFid_ = f->fid();
}


Deriv DerivOfSymbFunc::representMeAsFunctionalDeriv() const 
{
  const SymbolicFuncElement* f 
    = dynamic_cast<const SymbolicFuncElement*>(evaluatableArg());
  TEUCHOS_TEST_FOR_EXCEPTION(f==0, std::logic_error, "DerivOfSymbFunc::"
                     "representMeAsFunctionalDeriv(), 'this' pointer "
                     "is not a symbolic function");
  return Deriv(f, SpatialDerivSpecifier(mi()));
}

Evaluator* DerivOfSymbFunc::createEvaluator(const EvaluatableExpr* expr,
                                   const EvalContext& context) const
{
  return new DerivOfSymbFuncEvaluator(dynamic_cast<const DerivOfSymbFunc*>(expr), context);
}


bool DerivOfSymbFunc::lessThan(const ScalarExpr* other) const
{
  const DerivOfSymbFunc* d = dynamic_cast<const DerivOfSymbFunc*>(other);
  TEUCHOS_TEST_FOR_EXCEPTION(d==0, std::logic_error, "cast should never fail at this point");
  
  if (argFid_ < d->argFid_) return true;
  if (d->argFid_ < argFid_) return false;
  
  return DiffOp::lessThan(other);
}
