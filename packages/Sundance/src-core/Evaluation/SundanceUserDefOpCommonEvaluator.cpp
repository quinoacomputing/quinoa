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

#include "SundanceUserDefOpCommonEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceUserDefOp.hpp"
#include "SundanceUserDefOpElement.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


UserDefOpCommonEvaluator
::UserDefOpCommonEvaluator(const UserDefFunctor* functor,
                           const UserDefOpElement* expr,
                           const EvalContext& context)
  :  maxOrder_(0),
     argValueIndex_(functor->domainDim(), -1),
     argValueIsConstant_(functor->domainDim()),
     constArgDerivCache_(functor->rangeDim()),
     varArgDerivCache_(functor->rangeDim()),
     cacheIsValid_(false),
     functor_(functor)
{
  /* Find the indices to the zeroth derivative of each argument */
  for (int i=0; i<functor->domainDim(); i++)
    {
      const SparsitySuperset* sArg = expr->evaluatableChild(i)->sparsitySuperset(context).get();
      int numConst=0;
      int numVec=0;
      for (int j=0; j<sArg->numDerivs(); j++)
        {
          if (sArg->deriv(j).order() == 0) 
            {
              if (sArg->state(j)==VectorDeriv)
                {
                  argValueIndex_[i] = numVec;              
                }
              else
                {
                  argValueIndex_[i] = numConst;              
                }
              break;
            }
          if (sArg->state(j) == VectorDeriv) 
            {
              numVec++;
            }
          else
            {
              numConst++;
            }
        }
      /* Check to make sure a zeroth derivative has been found. */
      TEUCHOS_TEST_FOR_EXCEPTION(argValueIndex_[i]==-1, std::runtime_error,
                         "no zeroth derivative found for argument #" << i
                         << " of " << expr->toString());
    }
}




void UserDefOpCommonEvaluator
::evalAllComponents(const EvalManager& mgr,
                    const Array<RCP<Array<double> > >& constArgVals,
                    const Array<RCP<Array<RCP<EvalVector> > > >& vArgVals) const 
{
  Tabs tab0;
  int numPoints = EvalManager::stack().vecSize();
  SUNDANCE_MSG3(mgr.verb(), tab0 << "UDOpCommonEval::evalAllComponents()");
  SUNDANCE_MSG3(mgr.verb(), tab0 << "num points = " << numPoints);
  SUNDANCE_MSG2(mgr.verb(), tab0 << "max diff order = " << maxOrder_);

  TEUCHOS_TEST_FOR_EXCEPTION(numPoints==0, std::logic_error,
                     "Empty vector detected in evalArgDerivs()"); 

  /* Get an array of pointers for the argument vectors.
   * If some of the arguments are constant, copy them into vectors. */
  Array<RCP<EvalVector> > argVals(argValueIndex_.size());
  Array<double*> argPtrs(argValueIndex_.size());
  for (int q=0; q<argValueIndex_.size(); q++)
    {
      Tabs tab1;
      if (argValueIsConstant_[q]) 
        {
          argVals[q] = mgr.popVector();
          double* ptr = argVals[q]->start();
          double c =  (*(constArgVals[q]))[argValueIndex_[q]];
          for (int p=0; p<numPoints; p++)
            {
              ptr[p] = c;
            }
          argPtrs[q] = ptr;
        }
      else
        {
          argVals[q] = (*(vArgVals[q]))[argValueIndex_[q]];
          argPtrs[q] = argVals[q]->start();
        }
      SUNDANCE_MSG3(mgr.verb(), tab1 << "argument #" << q << " is:");
      Tabs tab2;
      SUNDANCE_MSG3(mgr.verb(), tab2 << argVals[q]->str());
    }

  /* Allocate vectors for the function values and derivatives */
  TEUCHOS_TEST_FOR_EXCEPTION(maxOrder_ > 2, std::runtime_error,
                     "Differentiation order " << maxOrder_ << ">2 not supported "
                     "for user-defined operators");
  int rangeDim = functor_->rangeDim();
  int domainDim = functor_->domainDim();
  int nTotal = 1;
  int numFirst = domainDim;
  int numSecond = domainDim*(domainDim+1)/2;
  if (maxOrder_ > 0) nTotal += numFirst;
  if (maxOrder_ > 1) nTotal += numSecond;
  int numResultVecs = nTotal * rangeDim;
  
  /* The resultVecs array contains pointers to the numerical vectors in the
   * cache of vector-valued arg derivs.
   */
  Array<double*> resultVecs(numResultVecs);

  
  for (int i=0; i<rangeDim; i++)
    {
      varArgDerivCache_[i].resize(nTotal);
      varArgDerivCache_[i][0] = mgr.popVector();
      varArgDerivCache_[i][0]->resize(numPoints);
      varArgDerivCache_[i][0]->setString(functor_->name(i));
      int d0Pos = i;
      SUNDANCE_MSG3(mgr.verb(), "zeroth deriv of elem #" << i << " is at " << d0Pos);
      resultVecs[d0Pos] = varArgDerivCache_[i][0]->start();
      if (maxOrder_ > 0)
        {
          int ptr = 0;
          for (int j=0; j<domainDim; j++)
            {
              varArgDerivCache_[i][j+1] = mgr.popVector();
              varArgDerivCache_[i][j+1]->resize(numPoints);
              int d1Pos = rangeDim + domainDim*i + j;
              SUNDANCE_MSG3(mgr.verb(), "first deriv (" << j << ") of elem #" << i 
                                 << " is at " << d1Pos);
              resultVecs[d1Pos] 
                = varArgDerivCache_[i][j+1]->start();
              varArgDerivCache_[i][j+1]->setString("D[" + functor_->name(i) 
                                                   + ", " 
                                                   + argVals[j]->str() + "]");
              if (maxOrder_ > 1)
                {
                  for (int k=0; k<=j; k++, ptr++)
                    {
                      int m = (1 + numFirst);
                      varArgDerivCache_[i][m+ptr] = mgr.popVector();
                      varArgDerivCache_[i][m+ptr]->resize(numPoints);
                      varArgDerivCache_[i][m+ptr]->setString("D[" 
                                                             + functor_->name(i) 
                                                             + ", {" 
                                                             + argVals[j]->str() 
                                                             + ", "
                                                             + argVals[k]->str() 
                                                             + "}]");
                      int d2Pos = rangeDim + domainDim*rangeDim 
                        + i*numSecond + ptr;
                      SUNDANCE_MSG3(mgr.verb(), "second deriv (" << j << ", " << k << 
                                         ") of elem #" << i << " is at " << d2Pos);
                      resultVecs[d2Pos] 
                        = varArgDerivCache_[i][m+ptr]->start();
                    }
                }
            }
        }
    }

  
  /* Call the user's callback function. The results will be written
   * into the cache of argument derivatives. */
  const double** in = const_cast<const double**>(&(argPtrs[0]));
  double** out = &(resultVecs[0]);

  functor_->evaluationCallback(numPoints, maxOrder_, in, out);

 

  markCacheAsValid();
}


