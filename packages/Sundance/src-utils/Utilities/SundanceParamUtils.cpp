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

#include "SundanceParamUtils.hpp"

using Teuchos::Array;
using Teuchos::ParameterList;

using std::ifstream;

namespace Sundance
{
ParameterList mergeParamLists(const ParameterList& pDef, 
  const ParameterList& pIn)
{
  ParameterList rtn = pDef;
  using namespace Teuchos;
  
  /* replace any defaults with overriden values */
  ParameterList::ConstIterator i;

  for (i=pDef.begin(); i!=pDef.end(); i++)
  {
    const ParameterEntry& eDef = pDef.entry(i);

    const std::string& name = pDef.name(i);
    const ParameterEntry* eIn = pIn.getEntryPtr(name);
    if (eIn != NULL)
    {
      if (eIn->isList() && eDef.isList())
      {
        ParameterList sub = mergeParamLists(
          getValue<ParameterList>(eDef),
          getValue<ParameterList>(*eIn));
        rtn.set(name, sub);

      }
      else if (eIn->isType<int>() && eDef.isType<int>())
      {
        rtn.set(name, Teuchos::any_cast<int>(eIn->getAny()));
      }
      else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(eIn->isList() && !eDef.isList(), 
          std::runtime_error, "mismatched parameters in mergeParams()");
        TEUCHOS_TEST_FOR_EXCEPTION(!eIn->isList() && eDef.isList(), 
          std::runtime_error, "mismatched parameters in mergeParams()");
        TEUCHOS_TEST_FOR_EXCEPT(1);
      }
    }
    else
    {
    }
  }
  return rtn;
}
}
