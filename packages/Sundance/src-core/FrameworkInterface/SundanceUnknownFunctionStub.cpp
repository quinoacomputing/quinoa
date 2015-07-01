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

#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



UnknownFunctionStub::UnknownFunctionStub(const std::string& name, 
  int tensorOrder,
  int dim, 
  const RCP<const UnknownFuncDataStub>& data)
  : SymbolicFunc(makeFuncID(tensorOrder), 
    rcp_dynamic_cast<const CommonFuncDataStub>(data)), data_(data)
{
  FunctionIdentifier myFid = fid();
  if (tensorOrder==0)
  {
    Expr u = new UnknownFuncElement(data, name, "", myFid);
    append(u);
  }
  else if (tensorOrder==1)
  {
    for (int d=0; d<dim; d++)
    {
      std::string suffix="[" + Teuchos::toString(d) + "]";
      FunctionIdentifier fid = myFid.createComponent(d);
      append(new UnknownFuncElement(data, name, suffix, fid));
    }
  }
  else 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "tensor order = " << tensorOrder
      << " not supported");
  }
}



UnknownFunctionStub::UnknownFunctionStub(const std::string& name, 
  const SpectralBasis& sbasis, int tensorOrder, int dim,
  const RCP<const UnknownFuncDataStub>& data)
  : SymbolicFunc(FunctionIdentifier(),
    rcp_dynamic_cast<const CommonFuncDataStub>(data)), data_(data)
{
  Array<FunctionIdentifier> cFid(sbasis.nterms());

  for (int n=0; n<sbasis.nterms(); n++)
  {
    cFid[n] = makeFuncID(tensorOrder);
  }
  
  if (tensorOrder==0 || dim==1)
  {
    Array<Expr> coeffs(sbasis.nterms());
    for (int n=0; n<sbasis.nterms(); n++)
    {
      std::string suffix="";
      if (sbasis.nterms()>1) suffix = "[" + Teuchos::toString(n) + "]";
      coeffs[n] = new UnknownFuncElement(data, name, suffix, cFid[n]);
    }
    append(new SpectralExpr(sbasis, coeffs));
  }
  else if (tensorOrder==1)
  {
    for (int d=0; d<dim; d++)
    {
      std::string suffix="[" + Teuchos::toString(d) + "]";
      Array<Expr> coeffs(sbasis.nterms());
      for (int n=0; n<sbasis.nterms(); n++)
      {
        FunctionIdentifier fid = cFid[n].createComponent(d);
        if (sbasis.nterms()>1) suffix += "[" + Teuchos::toString(n) + "]";
        coeffs[n]= new UnknownFuncElement(data, name, suffix, fid);
      }
      append(new SpectralExpr(sbasis, coeffs));
    }
  }
  else 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "tensor order = " << tensorOrder
      << " not supported");
  }
}


