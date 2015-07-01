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

#include "SundanceChainRuleSum.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvalVector.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceSet.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;


ChainRuleSum::ChainRuleSum(const MultipleDeriv& md,
                           int resultIndex,
                           bool resultIsConstant)
  : md_(md),
    resultIndex_(resultIndex),
    resultIsConstant_(resultIsConstant),
    argDerivIndex_(),
    argDerivIsConstant_(),
    terms_()
{;}


void ChainRuleSum::addTerm(int argDerivIndex, 
                           bool argDerivIsConstant,
                           const Array<DerivProduct>& sum)
{
  argDerivIndex_.append(argDerivIndex);
  argDerivIsConstant_.append(argDerivIsConstant);
  terms_.append(sum);
}


void ChainRuleSum
::evalConstant(const EvalManager& mgr,
               const Array<RCP<Array<double> > >& constantArgResults,
               const Array<double>& constantArgDerivs,
               double& result) const
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "ChainRuleSum::evalConstant()");
  result = 0.0;
  for (int i=0; i<numTerms(); i++)
    {
      const double& argDeriv = constantArgDerivs[argDerivIndex(i)];
      const Array<DerivProduct>& sumOfDerivProducts = terms(i);
      double innerSum = 0.0;
      for (int j=0; j<sumOfDerivProducts.size(); j++)
        {
          double prod = 1.0;
          const DerivProduct& p = sumOfDerivProducts[j];
          for (int k=0; k<p.numConstants(); k++)
            {
              const IndexPair& ip = p.constant(k);
              prod *= (*(constantArgResults[ip.argIndex()]))[ip.valueIndex()];
            }
          innerSum += prod;
        }
      result += innerSum*argDeriv;
    }
}


void ChainRuleSum
::evalVar(const EvalManager& mgr,
          const Array<RCP<Array<double> > >& constantArgResults,
          const Array<RCP<Array<RCP<EvalVector> > > > & vArgResults,
          const Array<double>& constantArgDerivs,
          const Array<RCP<EvalVector> >& varArgDerivs,
          RCP<EvalVector>& varResult) const
{
  Tabs tabs;
  SUNDANCE_VERB_HIGH(tabs << "ChainRuleSum::evalVar()");
  int vecSize=-1;
  for (int i=0; i<varArgDerivs.size(); i++)
    {
      int s = varArgDerivs[i]->length();
      TEUCHOS_TEST_FOR_EXCEPTION(vecSize != -1 && s != vecSize, std::logic_error,
                         "inconsistent vector sizes " << vecSize
                         << " and " << s);
      vecSize = s;
    } 
  for (int i=0; i<vArgResults.size(); i++)
    {
      for (int j=0; j<vArgResults[i]->size(); j++)
        {
          int s = (*(vArgResults[i]))[j]->length();
          TEUCHOS_TEST_FOR_EXCEPTION(vecSize != -1 && s != vecSize, std::logic_error,
                             "inconsistent vector sizes " << vecSize
                             << " and " << s);
          vecSize = s;
        }
    } 
  TEUCHOS_TEST_FOR_EXCEPT(vecSize==-1);
  
  varResult = mgr.popVector();
  varResult->resize(vecSize);
  varResult->setToConstant(0.0);

  for (int i=0; i<numTerms(); i++)
    {
      Tabs tab1;
      SUNDANCE_VERB_HIGH(tab1 << "term=" << i << " of " << numTerms());
      RCP<EvalVector> innerSum = mgr.popVector();
      innerSum->resize(vecSize);
      innerSum->setToConstant(0.0);
      const Array<DerivProduct>& sumOfDerivProducts = terms(i);

      SUNDANCE_VERB_HIGH(tab1 << "inner sum init = " << *innerSum
                         << ", num terms = " << terms(i).size());

      for (int j=0; j<sumOfDerivProducts.size(); j++)
        {
          Tabs tab2;
          SUNDANCE_VERB_HIGH(tab2 << "dp=" << j << " of " << sumOfDerivProducts.size());
          const DerivProduct& p = sumOfDerivProducts[j];
          double cc = p.coeff();
          SUNDANCE_VERB_HIGH(tab2 << "multiplicity=" << cc);
          for (int k=0; k<p.numConstants(); k++)
            {
              const IndexPair& ip = p.constant(k);
              cc *= (*(constantArgResults[ip.argIndex()]))[ip.valueIndex()];
            }
          if (p.numVariables()==0)
            {
              innerSum->add_S(cc);
            }
          else if (p.numVariables()==1)
            {
              const IndexPair& ip = p.variable(0);
              const EvalVector* v 
                = (*(vArgResults[ip.argIndex()]))[ip.valueIndex()].get();
              if (cc==1.0) innerSum->add_V(v);
              else innerSum->add_SV(cc, v);
            }
          else if (p.numVariables()==2)
            {
              const IndexPair& ip0 = p.variable(0);
              const EvalVector* v0 
                = (*(vArgResults[ip0.argIndex()]))[ip0.valueIndex()].get();
              const IndexPair& ip1 = p.variable(1);
              const EvalVector* v1
                = (*(vArgResults[ip1.argIndex()]))[ip1.valueIndex()].get();
              if (cc==1.0) innerSum->add_VV(v0, v1);
              else innerSum->add_SVV(cc, v0, v1);
            }
          else
            {
              const IndexPair& ip0 = p.variable(0);
              const EvalVector* v0 
                = (*(vArgResults[ip0.argIndex()]))[ip0.valueIndex()].get();
              RCP<EvalVector> tmp = v0->clone();
              for (int k=1; k<p.numVariables(); k++)
                {
                  const IndexPair& ip1 = p.variable(k);
                  const EvalVector* v1
                    = (*(vArgResults[ip1.argIndex()]))[ip1.valueIndex()].get();
                  tmp->multiply_V(v1);
                }
              if (cc==1.0) innerSum->add_V(tmp.get());
              else innerSum->add_SV(cc, tmp.get());
            }
          SUNDANCE_VERB_HIGH(tab2 << "inner sum=" << *innerSum);
        }

      int adi = argDerivIndex(i);
      if (argDerivIsConstant(i))
        {
          const double& df_dq = constantArgDerivs[adi];
          varResult->add_SV(df_dq, innerSum.get());
        }
      else
        {
          const EvalVector* df_dq = varArgDerivs[adi].get();
          SUNDANCE_VERB_HIGH(tab1 << "arg deriv=" << *df_dq);
          varResult->add_VV(df_dq, innerSum.get());
          SUNDANCE_VERB_HIGH(tab1 << "outer sum=" << *varResult);
        }
    }
}


