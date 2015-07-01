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

#include "SundanceSymbolicFunc.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceZeroExpr.hpp"

#include "SundanceDerivSet.hpp"
#include "PlayaTabs.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

SymbolicFunc::SymbolicFunc(const FunctionWithID& fid, 
    const RCP<const CommonFuncDataStub>& data)
  : ListExpr(), FunctionWithID(fid), commonData_(data)
{}


void SymbolicFunc::substituteZero() const 
{
  for (int i=0; i<this->size(); i++)
    {
      const SymbolicFuncElement* u 
        = dynamic_cast<const SymbolicFuncElement*>(element(i).ptr().get());
      TEUCHOS_TEST_FOR_EXCEPTION(u==0, std::logic_error, 
                         "Non-symbolic function "
                         << element(i).toString() 
                         << " detected in SymbolicFunc::substituteZero()");
      u->substituteZero();
    }
}

void SymbolicFunc
::substituteFunction(const RCP<DiscreteFunctionStub>& u0) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(this->size() != u0->size(), std::logic_error,
                     "Mismatch between sizes of symbolic " << toString()
                     << " and discrete func " << u0->toString()
                     << " in substituteFunction()");

  for (int i=0; i<this->size(); i++)
    {
      const SymbolicFuncElement* u 
        = dynamic_cast<const SymbolicFuncElement*>(element(i).ptr().get());
      TEUCHOS_TEST_FOR_EXCEPTION(u==0, std::logic_error, 
                         "Non-symbolic function "
                         << element(i).toString() 
                         << " detected in SymbolicFunc::substituteFunction()");

      RCP<DiscreteFuncElement> df 
        = rcp_dynamic_cast<DiscreteFuncElement>(u0->element(i).ptr());
      u->substituteFunction(df);
    }
}

