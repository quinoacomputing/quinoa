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

#include "SundanceListExpr.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

static Time& appendToListTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("append to list"); 
  return *rtn;
}

ListExpr::ListExpr()
  : ExprBase(), elements_()
{;}

ListExpr::ListExpr(const Array<Expr>& elements)
  : ExprBase(), elements_(elements)
{;}

void ListExpr::append(const Expr& expr)
{
  TimeMonitor timer(appendToListTimer());
  elements_.append(expr);
}

Expr ListExpr::flatten() const 
{
  Expr rtn = new ListExpr();

  for (int i=0; i<this->size(); i++)
    {
      Expr e = element(i).flatten();
      for (int j=0; j<e.size(); j++)
        {
          rtn.append(e[j]);
        }
    }

  return rtn;
}

Expr ListExpr::join(const Expr& other) const 
{
  Expr rtn = new ListExpr(elements_);
  
  for (int i=0; i<other.size(); i++)
    {
      rtn.append(other[i]);
    }

  return rtn;
}

int ListExpr::size() const
{
  return elements_.size();
}

int ListExpr::totalSize() const 
{
  int rtn = 0;

  for (int i=0; i<this->size(); i++)
    {
      rtn += elements_[i].totalSize();
    }

  return rtn;
}

std::ostream& ListExpr::toText(std::ostream& os, bool paren) const
{
  os << "{";
  for (int i=0; i<elements_.size(); i++)
    {
      elements_[i].ptr()->toText(os, paren);
      if (i < elements_.size()-1) os << ", ";
    }
  os << "}";
  return os;
}


XMLObject ListExpr::toXML() const 
{
  XMLObject rtn("ListExpr");
  for (int i=0; i<elements_.length(); i++)
    {
      rtn.addChild(elements_[i].toXML());
    }
  return rtn;
}


