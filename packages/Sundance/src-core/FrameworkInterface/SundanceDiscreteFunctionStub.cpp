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

#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"


using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



DiscreteFunctionStub::DiscreteFunctionStub(const std::string& name, 
  int tensorOrder,
  int dim, 
  const RCP<DiscreteFuncDataStub>& data,
  int listIndex)
  : ListExpr(), data_(data)
{
  initTensor(name, tensorOrder, dim, data,  listIndex);
}

void DiscreteFunctionStub::initTensor(const std::string& name, 
  int tensorOrder,
  int dim, 
  const RCP<DiscreteFuncDataStub>& data,
  int listIndex)
{
  FunctionIdentifier myFid = makeFuncID(tensorOrder);
  if (tensorOrder==0)
  {
    append(new DiscreteFuncElement(data, name, "", myFid, listIndex));
  }
  else if (tensorOrder==1)
  {
    for (int d=0; d<dim; d++)
    {
      std::string suffix="[" + Teuchos::toString(d) + "]";
      FunctionIdentifier fid = myFid.createComponent(d);
      append(new DiscreteFuncElement(data, name, suffix, fid, listIndex));
    }
  }
  else 
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "tensor order = " << tensorOrder
      << " not supported");
  }
}


void DiscreteFunctionStub::initTensorSpectral(const std::string& name, 
  const SpectralBasis& sbasis, 
  int tensorOrder,
  int dim, 
  const RCP<DiscreteFuncDataStub>& data,
  int listIndexOffset)
{
  Array<FunctionIdentifier> cFid(sbasis.nterms());
  Array<int> listIndex(sbasis.nterms());

  for (int n=0; n<sbasis.nterms(); n++)
  {
    cFid[n] = makeFuncID(tensorOrder);
    listIndex[n] = listIndexOffset + n;
  }
  
  if (tensorOrder==0 || dim==1)
  {
    Array<Expr> coeffs(sbasis.nterms());
    for (int n=0; n<sbasis.nterms(); n++)
    {
      std::string suffix="";
      if (sbasis.nterms()>1) suffix = "[" + Teuchos::toString(n) + "]";
      coeffs[n] = new DiscreteFuncElement(data, name, suffix, cFid[n], listIndex[n]);
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
        coeffs[n]= new DiscreteFuncElement(data, name, suffix, fid, listIndex[n]);
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




DiscreteFunctionStub::DiscreteFunctionStub(const Array<string>& name, 
  const Array<std::pair<int,int> >& tensorStructure,
  const RCP<DiscreteFuncDataStub>& data)
  : ListExpr(), data_(data)
{
  TEUCHOS_TEST_FOR_EXCEPT(name.size() != tensorStructure.size() && name.size()!=1);

  if (tensorStructure.size()==1)
  {
    int tensorOrder = tensorStructure[0].first;
    int dim = tensorStructure[0].second;
    initTensor(name[0], tensorOrder, dim, data, 0);
  }
  else
  {
    for (int i=0; i<tensorStructure.size(); i++)
    {
      std::string nm;
      if (name.size()==1) nm = name[0] + "[" + Teuchos::toString(i) + "]";
      else nm = name[i];
      append(new DiscreteFunctionStub(
               nm,
               tensorStructure[i].first,
               tensorStructure[i].second,
               data, i));
    }
  }
}



DiscreteFunctionStub::DiscreteFunctionStub(const std::string& name, 
  const SpectralBasis& sbasis, int tensorOrder, int dim,
  const RCP<DiscreteFuncDataStub>& data,
  int listIndex)
  : ListExpr(), data_(data)
{
  initTensorSpectral(name, sbasis, tensorOrder, dim, data, listIndex);
}



DiscreteFunctionStub::DiscreteFunctionStub(const Array<string>& name, 
  const SpectralBasis& sbasis,  
  const Array<std::pair<int,int> >& tensorStructure,
  const RCP<DiscreteFuncDataStub>& data)
  : ListExpr(), data_(data)
{
  TEUCHOS_TEST_FOR_EXCEPT(name.size() != tensorStructure.size());
   if (tensorStructure.size()==1)
  {
    int tensorOrder = tensorStructure[0].first;
    int dim = tensorStructure[0].second;
    initTensorSpectral(name[0], sbasis, tensorOrder, dim, data, 0);
  }
  else
  {
    for (int i=0; i<tensorStructure.size(); i++)
    {
      std::string nm;
      if (name.size()==1) nm = name[0] + "[" + Teuchos::toString(i) + "]";
      else nm = name[i];
      append(new DiscreteFunctionStub(
               nm,
               sbasis,
               tensorStructure[i].first,
               tensorStructure[i].second,
               data, i*sbasis.nterms()));
    }
  }
}


