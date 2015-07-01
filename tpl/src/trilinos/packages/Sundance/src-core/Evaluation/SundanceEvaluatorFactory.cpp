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

#include "SundanceEvaluatorFactory.hpp"
#include "SundanceInstructionCachingEvaluator.hpp"
#include "SundanceEvaluator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCurveNormExpr.hpp"
#include "SundanceCurveNormEvaluator.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"

using namespace Sundance;
using namespace Teuchos;


EvaluatorFactory::EvaluatorFactory()
{;}

Evaluator* EvaluatorFactory::commonCreate(const EvaluatableExpr* expr,
                                          const EvalContext& context,
                                          int topLevelDiffOrder) const
{
  const CoordExpr* c = dynamic_cast<const CoordExpr*>(expr);
  if (c != 0)
    {
      return new CoordExprEvaluator(c, context, topLevelDiffOrder);
    }
  
  const SpatiallyConstantExpr* sc 
    = dynamic_cast<const SpatiallyConstantExpr*>(expr);
  if (sc != 0)
    {
      return new ConstantEvaluator(sc, context, topLevelDiffOrder);
    }

  const SymbolicFuncElement* u
    = dynamic_cast<const SymbolicFuncElement*>(expr);
  if (u != 0)
    {
      return new SymbolicFuncElementEvaluator(u, context, topLevelDiffOrder);
    }
  
  const DiscreteFuncElement* df
    = dynamic_cast<const DiscreteFuncElement*>(expr);
  if (df != 0)
    {
      return new DiscreteFuncElementEvaluator(df, context, topLevelDiffOrder);
    }

  const CurveNormExpr* cne
    = dynamic_cast<const DiscreteFuncElement*>(expr);
  if (cne != 0)
    {
      return new CurveNormEvaluator(cne, context, topLevelDiffOrder);
    }

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                     "EvaluatorFactory::commonCreate() could not create an "
                     "evaluator for " << expr->toString());

  return 0;
}


RCP<EvaluatorFactory>&  EvaluatorFactory::defaultEvaluator()
{
  static RCP<EvaluatorFactory> rtn 
    = rcp(new InstructionCachingEvaluatorFactory());
  return rtn;
}
