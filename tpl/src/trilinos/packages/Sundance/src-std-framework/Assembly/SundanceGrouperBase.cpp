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

#include "SundanceGrouperBase.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceFuncWithBasis.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceUnknownParameterElement.hpp"
#include "SundanceTestFunction.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;



void GrouperBase::setVerb(
  int setupVerb,
  int integrationVerb,
  int transformVerb)
{
  setupVerb_ = setupVerb;
  integrationVerb_ = integrationVerb;
  transformVerb_ = transformVerb;
}



void GrouperBase::extractWeakForm(const EquationSet& eqn,
  const MultipleDeriv& functionalDeriv,
  BasisFamily& varBasis, 
  BasisFamily& unkBasis,
  MultiIndex& miVar, MultiIndex& miUnk,
  int& rawVarID, int& rawUnkID,  
  int& reducedVarID, int& reducedUnkID,  
  int& testBlock, int& unkBlock, 
  int& rawParamID, int& reducedParamID, 
  bool& isOneForm, bool& hasParam) const
{
  Tabs tab0(0);

  MultipleDeriv::const_iterator iter;

  isOneForm = false;  
  hasParam = false;

  if (functionalDeriv.size()==0) return;

  TEUCHOS_TEST_FOR_EXCEPTION(functionalDeriv.size() > 2, std::logic_error,
    "GrouperBase::extractWeakForm detected a functional "
    "derivative of order > 2: " 
    << functionalDeriv.toString());

  bool foundUnk = false;
  bool foundVar = false;

  SUNDANCE_MSG2(setupVerb(), 
    tab0 << "extracting weak form for functional derivative " 
    << functionalDeriv);


  for (iter = functionalDeriv.begin(); iter != functionalDeriv.end(); iter++)
  {
    Tabs tab;
    const Deriv& d = *iter;
      
    TEUCHOS_TEST_FOR_EXCEPTION(!d.isFunctionalDeriv(), std::logic_error,
      "GrouperBase::extractWeakForm "
      "detected a non-functional derivative: "
      << functionalDeriv.toString());
      
    const FunctionIdentifier& fid = d.fid();
      
    const SymbolicFuncElement* s = d.symbFuncElem();

    TEUCHOS_TEST_FOR_EXCEPTION(s==0, std::logic_error, 
      "GrouperBase::extractWeakForm failed to cast "
      "function to SymbolicFuncElement");
      

    int dofID = fid.dofID();
    int myIndex = fid.componentIndex();

    if (!foundVar && eqn.hasVarID(dofID))
    {
      TEUCHOS_TEST_FOR_EXCEPTION(d.isParameter(), std::logic_error,
        "Parameter not expected here");
      foundVar = true;
      reducedVarID = eqn.reducedVarID(dofID);
      rawVarID = dofID;
      testBlock = eqn.blockForVarID(dofID);

      SUNDANCE_MSG2(setupVerb(), 
        tab << "found varID=" << reducedVarID);

      const UnknownFuncElement* u
        = dynamic_cast<const UnknownFuncElement*>(s);

      const TestFuncElement* t
        = dynamic_cast<const TestFuncElement*>(s);

      TEUCHOS_TEST_FOR_EXCEPTION(u==0 && t==0, std::logic_error, 
        "GrouperBase::extractWeakForm could not cast "
        "variational function to either an "
        "UnknownFuncElement or a TestFuncElement");

      if (t != 0) 
      {
        varBasis = TestFunctionData::getData(t)->basis()[myIndex];
      }
      else 
      {
        varBasis = UnknownFunctionData::getData(u)->basis()[myIndex];
      }
      SUNDANCE_MSG2(setupVerb(), 
        tab << "found varBasis=" << varBasis);

      miVar = d.opOnFunc().mi();
      SUNDANCE_MSG2(setupVerb(), 
        tab << "found var multi index=" << miVar.toString());
    }
    else if (eqn.hasFixedParamID(dofID))
    {
      const UnknownParameterElement* upe
        = dynamic_cast<const UnknownParameterElement*>(s);
      TEUCHOS_TEST_FOR_EXCEPTION(upe==0, std::logic_error, 
        "GrouperBase::extractWeakForm could not cast "
        "unknown parameter to UnknownParameterElement");
      hasParam = true;
      rawParamID = dofID;
      reducedParamID = eqn.reducedFixedParamID(dofID);
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(d.isParameter(), std::logic_error,
        "Parameter not expected here");
      const UnknownFuncElement* u
        = dynamic_cast<const UnknownFuncElement*>(s);
      TEUCHOS_TEST_FOR_EXCEPTION(u==0, std::logic_error, 
        "GrouperBase::extractWeakForm could not cast "
        "unknown function to UnknownFuncElement");
      foundUnk = true;
      reducedUnkID = eqn.reducedUnkID(dofID);
      rawUnkID = dofID;
      unkBlock = eqn.blockForUnkID(dofID);

      SUNDANCE_MSG2(setupVerb(), 
        tab << "found reducedUnkID=" << reducedUnkID);

      unkBasis = UnknownFunctionData::getData(u)->basis()[myIndex];
      SUNDANCE_MSG2(setupVerb(), 
        tab << "found unkBasis=" << unkBasis);

      miUnk = d.opOnFunc().mi();
      SUNDANCE_MSG2(setupVerb(), 
        tab << "found unk multi index=" << miUnk.toString());
    }
  }

  if (!foundUnk) isOneForm = true;
}
