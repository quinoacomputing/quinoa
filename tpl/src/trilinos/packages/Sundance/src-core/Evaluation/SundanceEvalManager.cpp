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


#include "SundanceEvalManager.hpp"
#include "SundanceOut.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceCurveNormExpr.hpp"
#include "SundanceCellVectorExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceMultiIndex.hpp"



using namespace Sundance;
using namespace Teuchos;

EvalManager::EvalManager()
  : verb_(0),
    region_(),
    mediator_()
{}

void EvalManager::setVerb(int verb)
{
  verb_ = verb;
//  if (mediator_.get()) mediator_->setVerb(verb);
}

void EvalManager::evalCoordExpr(const CoordExpr* expr,
                                RCP<EvalVector>& result) const 
{

  TimeMonitor timer(coordEvalTimer());
  TEUCHOS_TEST_FOR_EXCEPTION(mediator() == 0, std::logic_error,
                     "uninitialized mediator in "
                     "EvalManager::evalCoordExpr");
  mediator()->evalCoordExpr(expr, result);
}


void EvalManager::evalCellDiameterExpr(const CellDiameterExpr* expr,
                                RCP<EvalVector>& result) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(mediator() == 0, std::logic_error,
                     "uninitialized mediator in "
                     "EvalManager::evalCellDiameterExpr");
  mediator()->evalCellDiameterExpr(expr, result);
}

void EvalManager::evalCurveNormExpr(const CurveNormExpr* expr,
                                RCP<EvalVector>& result) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(mediator() == 0, std::logic_error,
                     "uninitialized mediator in "
                     "EvalManager::evalCurveNormExpr");
  mediator()->evalCurveNormExpr(expr, result);
}

void EvalManager::evalCellVectorExpr(const CellVectorExpr* expr,
                                RCP<EvalVector>& result) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(mediator() == 0, std::logic_error,
                     "uninitialized mediator in "
                     "EvalManager::evalCellVectorExpr");
  mediator()->evalCellVectorExpr(expr, result);
}


void EvalManager::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                                          const Array<MultiIndex>& mi,
                                          Array<RCP<EvalVector> >& result) const 
{
  TimeMonitor timer(discFuncEvalTimer());

  TEUCHOS_TEST_FOR_EXCEPTION(mediator() == 0, std::logic_error,
                     "uninitialized mediator in "
                     "EvalManager::evalDiscreteFuncElement");

  mediator()->evalDiscreteFuncElement(expr, mi, result);
}

void EvalManager::showResults(std::ostream& os,
			      const RCP<SparsitySuperset>& sparsity,
			      const Array<RCP<EvalVector> >& vecResults,
			      const Array<double>& constantResults) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(mediator() == 0, std::logic_error,
                     "uninitialized mediator in "
                     "EvalManager::showResults");

  mediator()->showResults(os, sparsity, vecResults, constantResults);
}



RCP<EvalVector> EvalManager::popVector() const
{
  return stack().popVector();
}

TempStack& EvalManager::stack()
{
  static TempStack rtn;
  return rtn;
}
