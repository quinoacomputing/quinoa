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

#include "SundanceCellVectorExpr.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsitySuperset.hpp"
#include "SundanceOut.hpp"
#include "SundanceObjectWithVerbosity.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;



namespace Sundance
{

Expr CellNormalExpr(int dimension, const std::string& name)
{
  Array<Expr> comps(dimension);
  for (int i=0; i<dimension; i++)
  {
    comps[i] = new CellVectorExpr(i, dimension, name + "["
      + Teuchos::toString(i) + "]");
  }
  return new ListExpr(comps);
}


Expr CellTangentExpr(int dimension, const std::string& name)
{
  Array<Expr> comp(dimension);
  for (int i=0; i<dimension; i++)
  {
    comp[i] = new CellVectorExpr(0, i, dimension, name + "("
        + Teuchos::toString(i) + ")");
  }
  return new ListExpr(comp);
}

}

CellVectorExpr::CellVectorExpr(int tangentBasisIndex, 
			       int tangentComponentIndex,
			       int dim,
			       const std::string& name)
  : EvaluatableExpr(), name_(name), dim_(dim), type_(CellTangentSpace),
    basisMemberIndex_(tangentBasisIndex), 
    componentIndex_(tangentComponentIndex)
{}


CellVectorExpr::CellVectorExpr(int normalComponentIndex, int dim,
  const std::string& name)
  : EvaluatableExpr(), name_(name), dim_(dim), type_(CellNormalVector),
    basisMemberIndex_(-1), componentIndex_(normalComponentIndex)
{}

bool CellVectorExpr::lessThan(const ScalarExpr* other) const
{
  const CellVectorExpr* f = dynamic_cast<const CellVectorExpr*>(other);
  TEUCHOS_TEST_FOR_EXCEPTION(f==0, std::logic_error, "cast should never fail at this point");
  if (type_ < f->type_) return true;
  if (type_ > f->type_) return false;
  if (dim_ < f->dim_) return true;
  if (dim_ > f->dim_) return false;
  if (basisMemberIndex_ < f->basisMemberIndex_) return true;
  if (basisMemberIndex_ > f->basisMemberIndex_) return false;
  return componentIndex_ < f->componentIndex_;
}


XMLObject CellVectorExpr::toXML() const 
{
  XMLObject rtn("CellVectorExpr");
  rtn.addAttribute("name", name_);
  return rtn;
}



Set<MultipleDeriv> 
CellVectorExpr::internalFindW(int order, const EvalContext& context) const
{
  Set<MultipleDeriv> rtn;
  
  if (order==0) rtn.put(MultipleDeriv());
  
  return rtn;
}




std::ostream& CellVectorExpr::toText(std::ostream& os, bool paren) const
{
  os << name();
  return os;
}




