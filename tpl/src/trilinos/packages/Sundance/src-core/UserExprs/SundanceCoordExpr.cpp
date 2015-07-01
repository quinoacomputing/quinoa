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

#include "SundanceCoordExpr.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"
#include "SundanceObjectWithVerbosity.hpp"

using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


CoordExpr::CoordExpr(int dir, const std::string& name)
  : EvaluatableExpr(),
    dir_(dir),
    name_(coordName(dir, name))
{}

bool CoordExpr::lessThan(const ScalarExpr* other) const
{
  const CoordExpr* c = dynamic_cast<const CoordExpr*>(other);
  TEUCHOS_TEST_FOR_EXCEPTION(c==0, std::logic_error, "cast should never fail at this point");
  return dir() < c->dir();
}

XMLObject CoordExpr::toXML() const 
{
  XMLObject rtn("CoordExpr");
  rtn.addAttribute("dir", Teuchos::toString(dir_));
  rtn.addAttribute("name", name());
  return rtn;
}

string CoordExpr::coordName(int dir, const std::string& name)
{
  if (name.length() > 0) return name;
  switch(dir)
    {
    case 0:
      return "x";
    case 1:
      return "y";
    case 2:
      return "z";
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                         "CoordExpr::coordName direction out of range [0,2]");
      return "error";
    }
}


Set<MultipleDeriv> 
CoordExpr::internalFindW(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "CoordExpr::internalFindW() for " << toString());  
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());

  if (order==1) 
    {
      Deriv x = coordDeriv(dir_);
      MultipleDeriv md;
      md.put(x);
      rtn.put(md);
    }

  SUNDANCE_VERB_HIGH(tab0 << "W[" << order << "]=" << rtn);
  return rtn;
}


Set<MultipleDeriv> 
CoordExpr::internalFindV(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "CoordExpr::internalFindV() for " << toString());  
  Set<MultipleDeriv> rtn;

  if (order==0) rtn.put(MultipleDeriv());

  SUNDANCE_VERB_HIGH(tab0 << "V[" << order << "]=" << rtn);
  return rtn;
}


Set<MultipleDeriv> 
CoordExpr::internalFindC(int order, const EvalContext& context) const
{
  Tabs tab0;
  SUNDANCE_VERB_HIGH(tab0 << "CoordExpr::internalFindC() for " << toString());  
  Set<MultipleDeriv> rtn;

  if (order==1) 
    {
      Deriv x = coordDeriv(dir_);
      MultipleDeriv md;
      md.put(x);
      rtn.put(md);
    }
  SUNDANCE_VERB_HIGH(tab0 << "C[" << order << "]=" << rtn);
  return rtn;
}






