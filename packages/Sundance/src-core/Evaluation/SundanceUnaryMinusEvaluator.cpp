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

#include "SundanceUnaryMinusEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceUnaryMinus.hpp"

#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





UnaryMinusEvaluator
::UnaryMinusEvaluator(const UnaryMinus* expr,
                      const EvalContext& context)
  : UnaryEvaluator<UnaryMinus>(expr, context)
{
  try
    {
      int vecResultIndex = 0;
      int constResultIndex = 0;
      
      for (int i=0; i<this->sparsity()->numDerivs(); i++)
        {
          /* Determine the index into which the result will be written */
          bool resultIsConstant = this->sparsity()->state(i)==ConstantDeriv; 
          
          if (!resultIsConstant)
            {
              addVectorIndex(i, vecResultIndex);
              vecResultIndex++;
            }
          else
            {
              addConstantIndex(i, constResultIndex);
              constResultIndex++;
            }
        }
    }
  catch(std::exception& e)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
                         "exception detected in UnaryMinusEvaluator: expr="
                         << expr->toString() << std::endl
                         << "exception=" << e.what());
    }
}

void UnaryMinusEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RCP<EvalVector> >& vectorResults) const
{
  Tabs tab;
  SUNDANCE_MSG1(mgr.verb(),
    tab << "UnaryMinusEvaluator::eval() expr=" << expr()->toString());


  /* evaluate the argument */
  evalOperand(mgr, constantResults, vectorResults);


  if (mgr.verb() > 2)
    {
      Out::os() << tab << "UnaryMinus operand results" << std::endl;
      argSparsitySuperset()->print(Out::os(), vectorResults,
                           constantResults);
    }

  for (int i=0; i<constantResults.size(); i++)
    {
      constantResults[i] *= -1;
    }

  for (int i=0; i<vectorResults.size(); i++)
    {
      vectorResults[i]->multiply_S(-1.0);
    }

  
  if (mgr.verb() > 1)
    {
      Out::os() << tab << "UnaryMinus results" << std::endl;
      sparsity()->print(Out::os(), vectorResults,
                           constantResults);
    }

  
}


