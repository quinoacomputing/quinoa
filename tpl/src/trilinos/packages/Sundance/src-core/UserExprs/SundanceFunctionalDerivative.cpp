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

#include "SundanceFunctionalDerivative.hpp"

#include "SundanceSymbolicTransformation.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceExplicitFunctionalDerivativeElement.hpp"

using namespace Sundance;
using namespace Teuchos;

using namespace Sundance;

Expr Sundance::FunctionalDerivative(const Expr& F, const Expr& u)
{
  if (F.size() == 1 && u.size()==1)
  {
    RCP<ScalarExpr> arg = SymbolicTransformation::getScalar(F[0]);

    const UnknownFuncElement* uPtr
      = dynamic_cast<const UnknownFuncElement*>(u[0].ptr().get());
    
    Deriv fd = funcDeriv(uPtr);
    
    return new ExplicitFunctionalDerivativeElement(arg, fd);
  }
  else if (F.size() > 1)
  {
    Array<Expr> rtnList(F.size());
    for (int i=0; i<rtnList.size(); i++)
    {
      rtnList[i] = FunctionalDerivative(F[i], u);
    }
    return new ListExpr(rtnList);
  }
  else 
  {
    Array<Expr> rtnList(u.size());
    for (int i=0; i<rtnList.size(); i++)
    {
      rtnList[i] = FunctionalDerivative(F, u[i]);
    }
    return new ListExpr(rtnList);
  }
  
}

