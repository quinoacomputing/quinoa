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

#include "SundanceEFDEEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSet.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


EFDEEvaluator::EFDEEvaluator(
  const ExplicitFunctionalDerivativeElement* expr, 
  const EvalContext& context
  )
  : UnaryEvaluator<ExplicitFunctionalDerivativeElement>(expr, context),
    constValIndexToArgIndexMap_(),
    varValIndexToArgIndexMap_()
{

  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "initializing EFDE evaluator for " 
                    << expr->toString());
  SUNDANCE_VERB_MEDIUM(tabs << "return sparsity " << std::endl << *(this->sparsity)());

  /* 
   * This evaluator requires no calculations. All that is done is to
   * map derivatives (md, fd) in the argument's result arrays to 
   * derivatives (md) in this expression's result arrays. 
   */
  

  int vecResultIndex = 0;
  int constResultIndex = 0;

  const Deriv& fd = expr->fd();

  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = this->sparsity()->deriv(i);
      const DerivState& myState = this->sparsity()->state(i);

      if (myState==ConstantDeriv)
      {
        Tabs tab2;
        SUNDANCE_VERB_HIGH(tab2 
          << "deriv is constant, will be stored at index "
          << constResultIndex << " in the const result array");
        addConstantIndex(i, constResultIndex++);
      }
      else
      {
        Tabs tab2;
        SUNDANCE_VERB_HIGH(tab2 
          << "deriv is variable, will be stored at index "
          << vecResultIndex << " in the var result array");
        addVectorIndex(i, vecResultIndex++);
      }
      
      MultipleDeriv dArg = d;
      dArg.put(fd);

      int argIndex = argSparsitySuperset()->getIndex(dArg);

      
      TEUCHOS_TEST_FOR_EXCEPTION(argIndex==-1, std::runtime_error,
        "Derivative " << dArg << " expected in argument but not found");

      
      const DerivState& argState = argSparsitySuperset()->state(argIndex);
      TEUCHOS_TEST_FOR_EXCEPTION(argState != myState, std::logic_error, 
        "mismatched states");

      if (argState==ConstantDeriv)
      {
        int constArgIndex = argEval()->constantIndexMap().get(argIndex);
        constValIndexToArgIndexMap_.append(constArgIndex);
      }
      else
      {
        int vectorArgIndex = argEval()->vectorIndexMap().get(argIndex);
        varValIndexToArgIndexMap_.append(vectorArgIndex);
      }
    }
  
  SUNDANCE_VERB_HIGH(tabs 
    << " constant index map " 
    << constValIndexToArgIndexMap_ << std::endl 
    << " vector index map " 
    << varValIndexToArgIndexMap_
    );
}



void EFDEEvaluator::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const 
{
  TimeMonitor timer(efdeEvalTimer());
  Tabs tabs(0);

  SUNDANCE_MSG1(mgr.verb(), tabs << "EFDEEvaluator::eval() expr=" 
    << expr()->toString());

  SUNDANCE_MSG3(mgr.verb(), tabs << "sparsity = " << std::endl 
    << *(this->sparsity)());

  constantResults.resize(constValIndexToArgIndexMap_.size());
  vectorResults.resize(varValIndexToArgIndexMap_.size());

  /* evaluate the argument */
  Array<RCP<EvalVector> > argVectorResults;
  Array<double> argConstantResults;

  evalOperand(mgr, argConstantResults, argVectorResults);


  if (mgr.verb() > 2)
    {
      Tabs tab1;
      Out::os() << tab1 << "EFDE operand results" << std::endl;
      mgr.showResults(Out::os(), argSparsitySuperset(), argVectorResults,
		      argConstantResults);
    }


  for (int i=0; i<constantResults.size(); i++)
  {
    constantResults[i] = argConstantResults[constValIndexToArgIndexMap_[i]];
  }

  
  for (int i=0; i<vectorResults.size(); i++)
  {
    vectorResults[i] = mgr.popVector();
    const RCP<EvalVector>& v = argVectorResults[varValIndexToArgIndexMap_[i]];
    vectorResults[i]->setTo_V(v.get());
  }

  
  

  if (mgr.verb() > 2)
  {
    Tabs tab1;
    Out::os() << tab1 << "EFDE results " << std::endl;
    mgr.showResults(Out::os(), this->sparsity(), vectorResults,
      constantResults);
  }
}

