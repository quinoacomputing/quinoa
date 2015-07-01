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

#include "SundanceSumEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSumExpr.hpp"
#include "SundanceProductExpr.hpp"

#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;




SumEvaluator::SumEvaluator(const SumExpr* se,
                           const EvalContext& context)
  : BinaryEvaluator<SumExpr>(se, context), 
    sign_(se->sign()),
    singleRightConstant_(),
    singleRightVector_(),
    singleLeftConstant_(),
    singleLeftVector_(),
    ccSums_(),
    cvSums_(),
    vcSums_(),
    vvSums_()
{
  Tabs tabs(0);
  int verb = context.evalSetupVerbosity();

  if (verb > 1)
    {
      Out::os() << tabs << "initializing sum evaluator for " 
           << se->toString() << std::endl;
    }

  int constantCounter = 0;
  int vectorCounter = 0;

  if (verb > 2)
    {
      Out::os() << std::endl << tabs << "return sparsity ";
      this->sparsity()->print(Out::os());
      Out::os() << std::endl << tabs << "left sparsity ";
      leftSparsity()->print(Out::os());
      Out::os() << std::endl << tabs << "right sparsity " ;
      rightSparsity()->print(Out::os());
      Out::os() << std::endl;
      
      Out::os() << tabs << "left vector index map " << leftEval()->vectorIndexMap() << std::endl;
      Out::os() << tabs << "right vector index map " << rightEval()->vectorIndexMap() << std::endl;
      
      Out::os() << tabs << "left constant index map " << leftEval()->constantIndexMap() << std::endl;
      Out::os() << tabs << "right constant index map " << rightEval()->constantIndexMap() << std::endl;
    }

  for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      const MultipleDeriv& d = this->sparsity()->deriv(i);
      TEUCHOS_TEST_FOR_EXCEPTION(!leftSparsity()->containsDeriv(d) 
                         && !rightSparsity()->containsDeriv(d),
                         std::logic_error,
                         "deriv " << d.toString() 
                         << " was not found in either left or right operand "
                         "of expr " << expr()->toString());
      
      int iLeft = -1;
      int iRight = -1;
      
      if (leftSparsity()->containsDeriv(d)) 
        {
          iLeft = leftSparsity()->getIndex(d);
        }

      if (rightSparsity()->containsDeriv(d)) 
        {
          iRight = rightSparsity()->getIndex(d);
        }

      SUNDANCE_MSG2(verb, tabs << "deriv " << d);

      if (iLeft == -1) /* case where the left operand is zero */
        {
          SUNDANCE_MSG3(verb, tabs << "left operand is zero ");
          if (rightSparsity()->state(iRight)==ConstantDeriv)
            {
              int rc = rightEval()->constantIndexMap().get(iRight);
              singleRightConstant_.append(tuple(constantCounter, rc));
              addConstantIndex(i, constantCounter++);
              SUNDANCE_MSG4(verb, tabs << "single right constant " 
                                 << singleRightConstant_[singleRightConstant_.size()-1]);
            }
          else
            {
              int rv = rightEval()->vectorIndexMap().get(iRight);
              singleRightVector_.append(tuple(vectorCounter, rv));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_MSG4(verb, tabs << "single right vector " 
                                 << singleRightVector_[singleRightVector_.size()-1]);
            }
        }
      else if (iRight == -1) /* case where the right operand is zero */
        {
          SUNDANCE_MSG3(verb, tabs << "right operand is zero ");
          if (leftSparsity()->state(iLeft)==ConstantDeriv)
            {
              int lc = leftEval()->constantIndexMap().get(iLeft);
              singleLeftConstant_.append(tuple(constantCounter, lc));
              addConstantIndex(i, constantCounter++);
              SUNDANCE_MSG4(verb, tabs << "single left constant " 
                                 << singleLeftConstant_[singleLeftConstant_.size()-1]);
            }
          else
            {
              int lv = leftEval()->vectorIndexMap().get(iLeft);
              singleLeftVector_.append(tuple(vectorCounter, lv));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_MSG4(verb, tabs << "single left vector " 
                                 << singleLeftVector_[singleLeftVector_.size()-1]);
            }
        }
      else /* both are nonzero */
        {
          SUNDANCE_MSG4(verb, tabs << "both operands are nonzero ");
          bool leftIsConstant = leftSparsity()->state(iLeft)==ConstantDeriv;
          bool rightIsConstant = rightSparsity()->state(iRight)==ConstantDeriv;
          
          if (leftIsConstant && rightIsConstant)
            {
              SUNDANCE_MSG4(verb, tabs << "both operands are constant");
              int lc = leftEval()->constantIndexMap().get(iLeft);
              int rc = rightEval()->constantIndexMap().get(iRight);
              ccSums_.append(tuple(constantCounter, lc, rc));
              addConstantIndex(i, constantCounter++);
              SUNDANCE_MSG4(verb, tabs << "c-c sum " << ccSums_[ccSums_.size()-1]);
            }
          else if (leftIsConstant)
            {
              SUNDANCE_MSG4(verb, tabs << "left operand is constant");
              int lc = leftEval()->constantIndexMap().get(iLeft);
              int rv = rightEval()->vectorIndexMap().get(iRight);
              cvSums_.append(tuple(vectorCounter, lc, rv));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_MSG4(verb, tabs << "c-v sum " << cvSums_[cvSums_.size()-1]);
            }
          else if (rightIsConstant)
            {
              SUNDANCE_MSG4(verb, tabs << "right operand is constant");
              int lv = leftEval()->vectorIndexMap().get(iLeft);
              int rc = rightEval()->constantIndexMap().get(iRight);
              vcSums_.append(tuple(vectorCounter, lv, rc));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_MSG4(verb, tabs << "v-c sum " << vcSums_[vcSums_.size()-1]);
            }
          else
            {
              SUNDANCE_MSG4(verb, tabs << "both operands are vectors");
              int lv = leftEval()->vectorIndexMap().get(iLeft);
              int rv = rightEval()->vectorIndexMap().get(iRight);
              vvSums_.append(tuple(vectorCounter, lv, rv));
              addVectorIndex(i, vectorCounter++);
              SUNDANCE_MSG4(verb, tabs << "v-v sum " << vvSums_[vvSums_.size()-1]);
            }
        }
    }
}

void SumEvaluator
::internalEval(const EvalManager& mgr,
               Array<double>& constantResults,
               Array<RCP<EvalVector> >& vectorResults) const 
{ 
  //  TimeMonitor timer(evalTimer());
  Tabs tabs(0);

  SUNDANCE_MSG1(mgr.verb(),
               tabs << "SumEvaluator::eval() expr=" << expr()->toString());

  /* evaluate the children */
  Array<RCP<EvalVector> > leftVectorResults; 
  Array<RCP<EvalVector> > rightVectorResults; 
  Array<double> leftConstResults;
  Array<double> rightConstResults;
  evalChildren(mgr, leftConstResults, leftVectorResults,
               rightConstResults, rightVectorResults);

  if (verb() > 2)
    {
      Out::os() << tabs << "left operand " << std::endl;
      mgr.showResults(Out::os(), leftSparsity(), leftVectorResults,
                            leftConstResults);
      Out::os() << tabs << "right operand " << std::endl;
      mgr.showResults(Out::os(), rightSparsity(), rightVectorResults,
                             rightConstResults);
    }
  
  constantResults.resize(this->sparsity()->numConstantDerivs());
  vectorResults.resize(this->sparsity()->numVectorDerivs());

  /* Do constant terms with left=0 */
  for (int i=0; i<singleRightConstant_.size(); i++)
    {
      Tabs tab1;
      constantResults[singleRightConstant_[i][0]]
        = sign_*rightConstResults[singleRightConstant_[i][1]];
      SUNDANCE_MSG2(mgr.verb(), tab1 << "sum for "
                           << constantResultDeriv(singleRightConstant_[i][0])
                           << ": L=0, R=" 
                           << sign_*rightConstResults[singleRightConstant_[i][1]]);
    }

  /* Do constant terms with right=0 */
  for (int i=0; i<singleLeftConstant_.size(); i++)
    {
      Tabs tab1;
      constantResults[singleLeftConstant_[i][0]]
        = leftConstResults[singleLeftConstant_[i][1]];
      SUNDANCE_MSG2(mgr.verb(), tab1 << "sum for " 
                           << constantResultDeriv(singleLeftConstant_[i][0])
                           << ": L=" 
                           << leftConstResults[singleLeftConstant_[i][1]]
                           << " R=0");
    }

  /* Do vector terms with left=0 */
  for (int i=0; i<singleRightVector_.size(); i++)
    {
      Tabs tab1;
      if (sign_ < 0.0) rightVectorResults[singleRightVector_[i][1]]->multiply_S(sign_);
      vectorResults[singleRightVector_[i][0]]
        = rightVectorResults[singleRightVector_[i][1]];
      SUNDANCE_MSG2(mgr.verb(), tab1 << "sum for "
                           << vectorResultDeriv(singleRightVector_[i][0])
                           << ": L=0, R=" 
                           << sign_ << "*" << 
                           rightVectorResults[singleRightVector_[i][1]]->str());
    }

  /* Do vector terms with right=0 */
  for (int i=0; i<singleLeftVector_.size(); i++)
    { 
      Tabs tab1;
      vectorResults[singleLeftVector_[i][0]]
        = leftVectorResults[singleLeftVector_[i][1]];
      SUNDANCE_MSG2(mgr.verb(), tab1 << "sum for " 
                           << vectorResultDeriv(singleLeftVector_[i][0])
                           << ": L=" 
                           << leftVectorResults[singleLeftVector_[i][1]]->str()
                           << " R=0");
    }

  /** Do constant-constant terms */
  for (int i=0; i<ccSums_.size(); i++)
    {
      Tabs tab1;
      constantResults[ccSums_[i][0]]
        = leftConstResults[ccSums_[i][1]] 
        + sign_*rightConstResults[ccSums_[i][2]];
      SUNDANCE_MSG2(mgr.verb(), tab1 << "c-c sum for " 
                           << constantResultDeriv(ccSums_[i][0])
                           << ": L=" << leftConstResults[ccSums_[i][1]] 
                           << " R=" << sign_*rightConstResults[ccSums_[i][2]]);
    }

  /** Do constant-vector sums */
  for (int i=0; i<cvSums_.size(); i++)
    {
      Tabs tab1;
      RCP<EvalVector>& v = rightVectorResults[cvSums_[i][2]];
      SUNDANCE_MSG2(mgr.verb(), tab1 << "doing c-v sum for " 
                           << vectorResultDeriv(cvSums_[i][0])
                           << ": L=" << leftConstResults[cvSums_[i][1]] 
                           << " R=" << sign_ << "*" 
                           << rightVectorResults[cvSums_[i][2]]->str());
      if (isOne(sign_))
        {
          v->add_S(leftConstResults[cvSums_[i][1]]);
        }
      else
        {
          v->multiply_S_add_S(sign_, leftConstResults[cvSums_[i][1]]);
        }
      vectorResults[cvSums_[i][0]] = v;
    }

  /* Do vector-constant sums */
  for (int i=0; i<vcSums_.size(); i++)
    {
      Tabs tab1;
      RCP<EvalVector>& v = leftVectorResults[vcSums_[i][1]] ;
      SUNDANCE_MSG2(mgr.verb(), tab1 << "doing v-c sum for " 
                           << vectorResultDeriv(vcSums_[i][0])
                           << ": L=" << leftVectorResults[vcSums_[i][1]]->str()
                           << " R=" 
                           << sign_*rightConstResults[vcSums_[i][2]]);
      v->add_S(sign_*rightConstResults[vcSums_[i][2]]);
      vectorResults[vcSums_[i][0]] = v;
    }

  /* Do vector-vector sums */
  for (int i=0; i<vvSums_.size(); i++)
    {
      Tabs tab1;
      RCP<EvalVector>& v = leftVectorResults[vvSums_[i][1]];
      SUNDANCE_MSG2(mgr.verb(), tab1 << "doing v-v sum for " 
                           << vectorResultDeriv(vvSums_[i][0])
                           << ": L=" 
                           << leftVectorResults[vvSums_[i][1]]->str() 
                           << " R=" << sign_ << "*" 
                           << rightVectorResults[vvSums_[i][2]]->str());
      if (isOne(sign_))
        {
          v->add_V(rightVectorResults[vvSums_[i][2]].get());
        }
      else
        {
          v->add_SV(sign_, rightVectorResults[vvSums_[i][2]].get());
        }
      vectorResults[vvSums_[i][0]] = v;
    }
  
  if (mgr.verb() > 1)
    {
      Out::os() << tabs << "sum result " << std::endl;
      mgr.showResults(Out::os(), this->sparsity(), vectorResults,
		      constantResults);
    }
}


