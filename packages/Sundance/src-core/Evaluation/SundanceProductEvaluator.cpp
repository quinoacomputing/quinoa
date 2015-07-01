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

#include "SundanceProductEvaluator.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceProductExpr.hpp"

#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;





ProductEvaluator::ProductEvaluator(const ProductExpr* expr,
  const EvalContext& context)
  : BinaryEvaluator<ProductExpr>(expr, context),
    maxOrder_(this->sparsity()->maxOrder()),
    resultIndex_(maxOrder_+1),
    resultIsConstant_(maxOrder_+1),
    hasWorkspace_(maxOrder_+1),
    workspaceIsLeft_(maxOrder_+1),
    workspaceIndex_(maxOrder_+1),
    workspaceCoeffIndex_(maxOrder_+1),
    workspaceCoeffIsConstant_(maxOrder_+1),
    ccTerms_(maxOrder_+1),
    cvTerms_(maxOrder_+1),
    vcTerms_(maxOrder_+1),
    vvTerms_(maxOrder_+1),
    startingVectors_(maxOrder_+1),
    startingParities_(maxOrder_+1)
{
  int verb = context.evalSetupVerbosity();

  try
  {
    Tabs tabs(0);
    {
      Tabs tab;
      Tabs tabz;
      SUNDANCE_MSG1(verb,
        tabs << "initializing product evaluator for " 
        << expr->toString());

      SUNDANCE_MSG2(verb,
        tab << "return sparsity " << std::endl
        << tabz << *(this->sparsity)());

      SUNDANCE_MSG2(verb,
        tab << "left sparsity " << std::endl 
        << tabz << *(leftSparsity()) << std::endl
        << tabz << "right sparsity " << std::endl 
        << tabz << *(rightSparsity()));
  
      SUNDANCE_MSG3(verb,
        tab << "left vector index map " 
        << leftEval()->vectorIndexMap() << std::endl
        << tabz << "right vector index map " 
        << rightEval()->vectorIndexMap() << std::endl
        << tabz << "left constant index map " 
        << leftEval()->constantIndexMap() << std::endl
        << tabz << "right constant index map " 
        << rightEval()->constantIndexMap());
    }                
    int vecResultIndex = 0;
    int constResultIndex = 0;

    for (int i=0; i<this->sparsity()->numDerivs(); i++)
    {
      Tabs tab0;
      const MultipleDeriv& d = this->sparsity()->deriv(i);

      SUNDANCE_MSG2(verb,
        tabs << std::endl 
        << tabs << "finding rules for deriv " << d);

      int order = d.order();

      /* Determine the index into which the result will be written */
      bool resultIsConstant = this->sparsity()->state(i)==ConstantDeriv; 
      resultIsConstant_[order].append(resultIsConstant);
      if (resultIsConstant)
      {
        SUNDANCE_MSG3(verb,
          tab0 << std::endl 
          << tab0 << "result will be in constant index " << constResultIndex);
        resultIndex_[order].append(constResultIndex);
        addConstantIndex(i, constResultIndex);
        constResultIndex++;
      }
      else
      {
        SUNDANCE_MSG3(verb,
          tab0 << std::endl 
          << tab0 << "result will be in constant index " << vecResultIndex);
        resultIndex_[order].append(vecResultIndex);
        addVectorIndex(i, vecResultIndex);
        vecResultIndex++;
      }

      /* If possible, we want to do the calculations in-place, writing into
       * one of the two operand's results vectors for the same derivative.
       * Provided that we process derivatives in descending order, it is
       * safe to destroy the operands' result vectors.
       */
       
      int dnLeftIndex = leftSparsity()->getIndex(d);
      int dnRightIndex = rightSparsity()->getIndex(d);

      
      bool hasVectorWorkspace = false;
      bool workspaceIsLeft = false;
      int workspaceIndex = -1;
      int workspaceCoeffIndex = -1;
      bool workspaceCoeffIsConstant = false;


      bool dnLeftIsConst = false;
      if (dnLeftIndex != -1) 
      {
        dnLeftIsConst = leftSparsity()->state(dnLeftIndex)==ConstantDeriv;
      }
      bool dnRightIsConst = false;
      if (dnRightIndex != -1)
      {
        dnRightIsConst = rightSparsity()->state(dnRightIndex)==ConstantDeriv;
      }

      /* First try to use the left result as a workspace */
      if (dnLeftIndex != -1 && !dnLeftIsConst 
        && rightSparsity()->containsDeriv(MultipleDeriv()))
      {
        /* We can write onto left vector */
        hasVectorWorkspace = true;
        workspaceIndex = leftEval()->vectorIndexMap().get(dnLeftIndex);       
        SUNDANCE_MSG3(verb,
          tab0 << "using left as workspace");
        workspaceIsLeft = true;
        int d0RightIndex = rightSparsity()->getIndex(MultipleDeriv());
        bool d0RightIsConst = rightSparsity()->state(d0RightIndex)==ConstantDeriv;
        workspaceCoeffIsConstant = d0RightIsConst;
        if (d0RightIsConst)
        {
          workspaceCoeffIndex 
            = rightEval()->constantIndexMap().get(d0RightIndex);
        }
        else
        {
          workspaceCoeffIndex 
            = rightEval()->vectorIndexMap().get(d0RightIndex);
        }
      }

      /* If the left didn't provide a workspace, try the right */
      if (!hasVectorWorkspace && dnRightIndex != -1 && !dnRightIsConst
        && leftSparsity()->containsDeriv(MultipleDeriv()))
      {
        /* We can write onto right vector */
        hasVectorWorkspace = true;
        workspaceIndex = rightEval()->vectorIndexMap().get(dnRightIndex); 
        workspaceIsLeft = false;
        SUNDANCE_MSG3(verb,
          tab0 << "using right as workspace");
        int d0LeftIndex = leftSparsity()->getIndex(MultipleDeriv());
        bool d0LeftIsConst = leftSparsity()->state(d0LeftIndex)==ConstantDeriv;
        workspaceCoeffIsConstant = d0LeftIsConst;
        if (d0LeftIsConst)
        {
          workspaceCoeffIndex 
            = leftEval()->constantIndexMap().get(d0LeftIndex);
        }
        else
        {
          workspaceCoeffIndex 
            = leftEval()->vectorIndexMap().get(d0LeftIndex);
        }
      }
      
      if (hasVectorWorkspace && verb > 2)
      {
        std::string wSide = "right";
        MultipleDeriv wsDeriv;
        if (workspaceIsLeft) 
        {
          wSide = "left";
          wsDeriv = leftSparsity()->deriv(dnLeftIndex);
        }
        else
        {
          wsDeriv = rightSparsity()->deriv(dnRightIndex);
        }
        SUNDANCE_MSG2(verb, tab0 << "has workspace vector: "
          << wSide << " deriv= " 
          << wsDeriv
          << ", coeff index= " << workspaceCoeffIndex);
      }
      else
      {
        SUNDANCE_MSG2(verb, tab0 << "has no workspace vector");
      }

      hasWorkspace_[order].append(hasVectorWorkspace);
      workspaceIsLeft_[order].append(workspaceIsLeft);
      workspaceIndex_[order].append(workspaceIndex);
      workspaceCoeffIndex_[order].append(workspaceCoeffIndex);
      workspaceCoeffIsConstant_[order].append(workspaceCoeffIsConstant);
      
      ProductRulePerms perms;
      d.productRulePermutations(perms);
      SUNDANCE_MSG4(verb,
        tabs << "product rule permutations " << perms);

      Array<Array<int> > ccTerms;
      Array<Array<int> > cvTerms;
      Array<Array<int> > vcTerms;
      Array<Array<int> > vvTerms;
      Array<int> startingVector;
      ProductParity startingParity;

      bool hasStartingVector = false;

      for (ProductRulePerms::const_iterator 
             iter = perms.begin(); iter != perms.end(); iter++)
      {
        Tabs tab1;
        MultipleDeriv dLeft = iter->first.first();
        MultipleDeriv dRight = iter->first.second();
        int multiplicity = iter->second;
          
        /* skip this combination if we've already pulled it out 
         * for workspace */
        if (hasVectorWorkspace && workspaceIsLeft 
          && dLeft.order()==order) continue;
        if (hasVectorWorkspace && !workspaceIsLeft 
          && dRight.order()==order) continue;

        if (!leftSparsity()->containsDeriv(dLeft)
          || !rightSparsity()->containsDeriv(dRight)) continue;
        int iLeft = leftSparsity()->getIndex(dLeft);
        int iRight = rightSparsity()->getIndex(dRight);

        if (iLeft==-1 || iRight==-1) continue;
        SUNDANCE_MSG4(verb,
          tab1 << "left deriv=" << dLeft);
        SUNDANCE_MSG4(verb,
          tab1 << "right deriv=" << dRight);

        bool leftIsConst = leftSparsity()->state(iLeft)==ConstantDeriv;
        bool rightIsConst = rightSparsity()->state(iRight)==ConstantDeriv;

        if (leftIsConst)
        {
          Tabs tab2;
          SUNDANCE_MSG4(verb,
            tab2 << "left is const");
          int leftIndex = leftEval()->constantIndexMap().get(iLeft);
          if (rightIsConst)
          {
            SUNDANCE_MSG4(verb,
              tab2 << "right is const");
            int rightIndex = rightEval()->constantIndexMap().get(iRight);
            ccTerms.append(tuple(leftIndex, rightIndex, multiplicity));
          }
          else
          {
            SUNDANCE_MSG4(verb,
              tab2 << "right is vec");
            int rightIndex = rightEval()->vectorIndexMap().get(iRight);
            if (!hasVectorWorkspace && !hasStartingVector)
            {
              SUNDANCE_MSG4(verb,
                tab1 << "found c-v starting vec");
              startingVector = tuple(leftIndex, rightIndex, 
                multiplicity);
              hasStartingVector = true;
              startingParity = ConstVec;
            }
            else
            {
              SUNDANCE_MSG4(verb,
                tab1 << "found c-v term");
              cvTerms.append(tuple(leftIndex, rightIndex, 
                  multiplicity));
            }
          }
        }
        else
        {
          Tabs tab2;
          SUNDANCE_MSG4(verb,
            tab2 << "left is vec");
          int leftIndex = leftEval()->vectorIndexMap().get(iLeft);
          if (rightIsConst)
          {
            SUNDANCE_MSG4(verb,
              tab2 << "right is const");
            int rightIndex = rightEval()->constantIndexMap().get(iRight);
            if (!hasVectorWorkspace && !hasStartingVector)
            {
              SUNDANCE_MSG4(verb,
                tab1 << "found v-c starting vec");
              startingVector = tuple(leftIndex, rightIndex, 
                multiplicity);
              hasStartingVector = true;
              startingParity = VecConst;
            }
            else
            {
              SUNDANCE_MSG4(verb,
                tab1 << "found v-c term");
              vcTerms.append(tuple(leftIndex, rightIndex, 
                  multiplicity));
            }
          }
          else
          {
            SUNDANCE_MSG4(verb,
              tab2 << "right is vec");
            int rightIndex = rightEval()->vectorIndexMap().get(iRight);
            if (!hasVectorWorkspace && !hasStartingVector)
            {
              SUNDANCE_MSG4(verb,
                tab1 << "found v-v starting vec");
              startingVector = tuple(leftIndex, rightIndex, 
                multiplicity);
              hasStartingVector = true;
              startingParity = VecVec;
            }
            else
            {
              SUNDANCE_MSG4(verb,
                tab1 << "found v-v term");
              vvTerms.append(tuple(leftIndex, rightIndex, 
                  multiplicity));
            }
          }
        }
      }
      ccTerms_[order].append(ccTerms);
      cvTerms_[order].append(cvTerms);
      vcTerms_[order].append(vcTerms);
      vvTerms_[order].append(vvTerms);
      startingVectors_[order].append(startingVector);
      startingParities_[order].append(startingParity);

      if (verb > 2)
      {
        Tabs tab0;
        Out::os() << tab0 << "deriv " << i << " order=" << order ;
        if (resultIsConstant)
        {
          Out::os() << " constant result, index= ";
        }
        else
        {
          Out::os() << " vector result, index= ";
        }
        Out::os() << resultIndex_[order] << std::endl;
        {
          Tabs tab1;
          if (hasStartingVector) 
          {
            Out::os() << tab1 << "starting vector " << startingVector << std::endl;
          }
          Out::os() << tab1 << "c-c terms " << ccTerms << std::endl;
          Out::os() << tab1 << "c-v terms " << cvTerms << std::endl;
          Out::os() << tab1 << "v-c terms " << vcTerms << std::endl;
          Out::os() << tab1 << "v-v terms " << vvTerms << std::endl;
        }
      }
    }

    if (verb > 2)
    {
      Out::os() << tabs << "maps: " << std::endl;
      Out::os() << tabs << "vector index map " << vectorIndexMap() << std::endl;
      Out::os() << tabs << "constant index map " << constantIndexMap() << std::endl;
    }
  }
  catch(std::exception& e)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
      "exception detected in ProductEvaluator: expr="
      << expr->toString() << std::endl
      << "exception=" << e.what());
  }
}

void ProductEvaluator
::internalEval(const EvalManager& mgr,
  Array<double>& constantResults,
  Array<RCP<EvalVector> >& vectorResults) const
{
  Tabs tabs(0);

  SUNDANCE_MSG1(mgr.verb(),
    tabs << "ProductEvaluator::eval() expr=" 
    << expr()->toString());

  /* evaluate the children */
  Array<RCP<EvalVector> > leftVectorResults; 
  Array<RCP<EvalVector> > rightVectorResults; 
  Array<double> leftConstantResults;
  Array<double> rightConstantResults;
  evalChildren(mgr, leftConstantResults, leftVectorResults,
    rightConstantResults, rightVectorResults);

  if (mgr.verb() > 2)
  {
    Out::os() << tabs << "left operand results" << std::endl;
    mgr.showResults(Out::os(), leftSparsity(), leftVectorResults,
      leftConstantResults);
    Out::os() << tabs << "right operand results" << std::endl;
    mgr.showResults(Out::os(), rightSparsity(), rightVectorResults,
      rightConstantResults);
  }
  
  constantResults.resize(this->sparsity()->numConstantDerivs());
  vectorResults.resize(this->sparsity()->numVectorDerivs());

  /* Evaluate from high to low differentiation order. This is necessary
   * so that we can do in-place evaluation of derivatives, overwriting
   * one of the operands' vectors */

  for (int order=maxOrder_; order>=0; order--)
  {
    for (int i=0; i<resultIndex_[order].size(); i++)
    {
      double constantVal = 0.0;
      const Array<Array<int> >& ccTerms = ccTerms_[order][i];
      for (int j=0; j<ccTerms.size(); j++)
      {
        /* add L*R*multiplicity to result */
        constantVal += leftConstantResults[ccTerms[j][0]] 
          * rightConstantResults[ccTerms[j][1]] * ccTerms[j][2];
      }
      /* If this derivative is a constant, we're done. */
      if (resultIsConstant_[order][i]) 
      {
        constantResults[resultIndex_[order][i]] = constantVal;
        continue;
      }


      /* For the vector case, the first thing to do is get storage
       * for the result. If we can reuse one of the operands' results
       * as workspace, great. If not, we'll need to allocate a new
       * vector. */
      RCP<EvalVector> result;
      if (hasWorkspace_[order][i])
      {
        int index = workspaceIndex_[order][i];
        int coeffIndex = workspaceCoeffIndex_[order][i];
        bool coeffIsConstant = workspaceCoeffIsConstant_[order][i];
        if (workspaceIsLeft_[order][i])
        {
          result = leftVectorResults[index];
          if (coeffIsConstant)
          {
            const double& coeff = rightConstantResults[coeffIndex];
            if (!isOne(coeff)) result->multiply_S(coeff);
          }
          else
          {
            const RCP<EvalVector>& coeff 
              = rightVectorResults[coeffIndex];
            result->multiply_V(coeff.get());
          }
        }
        else
        {
          result = rightVectorResults[index];
          if (coeffIsConstant)
          {
            const double& coeff = leftConstantResults[coeffIndex];
            if (!isOne(coeff)) result->multiply_S(coeff);
          }
          else
          {
            const RCP<EvalVector>& coeff 
              = leftVectorResults[coeffIndex];
            result->multiply_V(coeff.get());
          }
        }
        SUNDANCE_MSG5(mgr.verb(), tabs << "workspace = " << result->str());
      }
      else
      {
        result = mgr.popVector();
        const Array<int>& sv = startingVectors_[order][i];
        int multiplicity = sv[2];
        switch(startingParities_[order][i])
        {
          case VecVec:
            if (!isZero(constantVal))
            {
              if (!isOne(multiplicity))
              {
                result->setTo_S_add_SVV(constantVal,
                  multiplicity,
                  leftVectorResults[sv[0]].get(),
                  rightVectorResults[sv[1]].get());
              }
              else
              {
                result->setTo_S_add_VV(constantVal,
                  leftVectorResults[sv[0]].get(),
                  rightVectorResults[sv[1]].get());
              }
            }
            else
            {
              if (!isOne(multiplicity))
              {
                result->setTo_SVV(multiplicity,
                  leftVectorResults[sv[0]].get(),
                  rightVectorResults[sv[1]].get());
              }
              else
              {
                result->setTo_VV(leftVectorResults[sv[0]].get(),
                  rightVectorResults[sv[1]].get());
              }
            }
            SUNDANCE_MSG5(mgr.verb(), tabs << "init to v-v prod, m="
              << multiplicity << ", left=" 
              << leftVectorResults[sv[0]]->str()
              << ", right=" 
              << rightVectorResults[sv[1]]->str()
              << ", result=" << result->str());
            break;
          case VecConst:
            if (!isZero(constantVal))
            {
              if (!isOne(multiplicity*rightConstantResults[sv[1]]))
              {
                result->setTo_S_add_SV(constantVal,
                  multiplicity
                  *rightConstantResults[sv[1]],  
                  leftVectorResults[sv[0]].get());
              }
              else
              {
                result->setTo_S_add_V(constantVal,
                  leftVectorResults[sv[0]].get());
              }
            }
            else
            {
              if (!isOne(multiplicity*rightConstantResults[sv[1]]))
              {
                result->setTo_SV(multiplicity
                  *rightConstantResults[sv[1]],  
                  leftVectorResults[sv[0]].get());
              }
              else
              {
                result->setTo_V(leftVectorResults[sv[0]].get());
              }
            }
            SUNDANCE_MSG5(mgr.verb(), tabs << "init to v-c prod, m="
              << multiplicity << ", left=" 
              << leftVectorResults[sv[0]]->str()
              << ", right=" 
              << rightConstantResults[sv[1]]
              << ", result=" << result->str());
            break;
          case ConstVec:
            if (!isZero(constantVal))
            {
              if (!isOne(multiplicity*leftConstantResults[sv[0]]))
              {
                result->setTo_S_add_SV(constantVal,
                  multiplicity
                  *leftConstantResults[sv[0]],  
                  rightVectorResults[sv[1]].get());
              }
              else
              {
                result->setTo_S_add_V(constantVal,  
                  rightVectorResults[sv[1]].get());
              }
            }
            else
            {
              if (!isOne(multiplicity*leftConstantResults[sv[0]]))
              {
                result->setTo_SV(multiplicity
                  *leftConstantResults[sv[0]],  
                  rightVectorResults[sv[1]].get());
              }
              else
              {
                result->setTo_V(rightVectorResults[sv[1]].get());
              }
            }
            SUNDANCE_MSG5(mgr.verb(), tabs << "init to c-v prod, m="
              << multiplicity << ", left=" 
              << leftConstantResults[sv[0]]
              << ", right=" 
              << rightVectorResults[sv[1]]->str()
              << ", result=" << result->str());
                  
        }
        SUNDANCE_MSG4(mgr.verb(), tabs << "starting value = " << result->str());
      }
      vectorResults[resultIndex_[order][i]] = result;

      const Array<Array<int> >& cvTerms = cvTerms_[order][i];
      const Array<Array<int> >& vcTerms = vcTerms_[order][i];
      const Array<Array<int> >& vvTerms = vvTerms_[order][i];

      for (int j=0; j<cvTerms.size(); j++)
      {
        SUNDANCE_MSG4(mgr.verb(), tabs << "adding c-v term " << cvTerms[j]);
              
        double multiplicity = cvTerms[j][2];
        double scalar = multiplicity*leftConstantResults[cvTerms[j][0]];
        const EvalVector* vec = rightVectorResults[cvTerms[j][1]].get();
        if (!isOne(scalar))
        {
          result->add_SV(scalar, vec);
        }
        else
        {
          result->add_V(vec);
        }
      } 

      for (int j=0; j<vcTerms.size(); j++)
      {
        SUNDANCE_MSG4(mgr.verb(), tabs << "adding v-c term " << vcTerms[j]);

        double multiplicity = vcTerms[j][2];
        double scalar = multiplicity*rightConstantResults[vcTerms[j][1]];
        const EvalVector* vec = leftVectorResults[vcTerms[j][0]].get();
        if (!isOne(scalar))
        {
          result->add_SV(scalar, vec);
        }
        else
        {
          result->add_V(vec);
        }
      }

      for (int j=0; j<vvTerms.size(); j++)
      {
        SUNDANCE_MSG4(mgr.verb(), tabs << "adding v-v term " << vvTerms[j]);

        double multiplicity = vvTerms[j][2];
        const EvalVector* vec1 = leftVectorResults[vvTerms[j][0]].get();
        const EvalVector* vec2 = rightVectorResults[vvTerms[j][1]].get();
        if (!isOne(multiplicity))
        {
          result->add_SVV(multiplicity, vec1, vec2);
        }
        else
        {
          result->add_VV(vec1, vec2);
        }
      }

          
    }
  }

  if (mgr.verb() > 1)
  {
    Out::os() << tabs << "product result " << std::endl;
    mgr.showResults(Out::os(), this->sparsity(), vectorResults,
      constantResults);
  }
}
