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

#include "SundanceFuncWithBasis.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"


namespace Sundance
{
using namespace Teuchos;

std::string describeFunction(const Expr& f)
{
  TEUCHOS_TEST_FOR_EXCEPT(f.ptr().get()==0);

  if (f.size() == 1)
  {
    const FuncElementBase* fe = dynamic_cast<const FuncElementBase*>(f[0].ptr().get());
    TEUCHOS_TEST_FOR_EXCEPTION(fe==0, std::runtime_error, "expected a FuncElementBase, "
      "found " << typeid(*fe).name());
    
    const UnknownFuncElement* u = dynamic_cast<const UnknownFuncElement*>(f[0].ptr().get());

    const TestFuncElement* t = dynamic_cast<const TestFuncElement*>(f[0].ptr().get());

    const DiscreteFuncElement* d = dynamic_cast<const DiscreteFuncElement*>(f[0].ptr().get());

    std::string type;
    if (t != 0) 
    {
      type = "TFElem";
    }
    else if (u != 0) 
    {
      type = "UFElem";
    }
    else if (d != 0)
    {
      type = "DFElem";
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "unrecognized function " 
        << f[0]);
    }

    std::string rtn = type + "[name=" + fe->name() + ", fid=" + fe->fid().toString() + "]";
    return rtn;
      
  }
  else
  {
    std::string rtn = "{";
    for (int i=0; i<f.size(); i++)
    {
      if (i != 0) rtn += ", ";
      rtn += describeFunction(f[i]);
    }
    rtn += "}";
    return rtn;
  }
}
}
