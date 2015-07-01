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


#include "SundanceRivaraEdge.hpp"
#include "SundanceRivaraElement.hpp"
#include "SundanceRivaraNode.hpp"
#include "SundanceRivaraMesh.hpp"

using namespace Sundance::Rivara;
using namespace Teuchos;



Edge::Edge(const RCP<Node>& a,
           const RCP<Node>& b)
  : label_(-1),nodes_(tuple(a,b)), elements_(), midpoint_(),
		ownerProc_()
{
	if (a->ownerProc() > b->ownerProc())
		{
			ownerProc_ = a->ownerProc();
		}
	else
		{
			ownerProc_ = b->ownerProc();
		}
}

void Edge::addConnectingElement(Element* tri)
{
  elements_.append(tri);
}

double Edge::length() const 
{
  const Point& x1 = nodes_[0]->pt();
  const Point& x2 = nodes_[1]->pt();

  return sqrt((x1-x2)*(x1-x2));
}

void Edge::getUnrefinedCofacets(Array<Element*>& c) const 
{
  for (int i=0; i<elements_.length(); i++)
    {
      if (!elements_[i]->hasChildren()) c.append(elements_[i]);
    }
}

RCP<Node> Edge::bisect(RivaraMesh* mesh)
{
  /* if we've already been bisected, return the existing midpoint node */
  if (!(midpoint_.get() == 0))
    {
      return midpoint_;
    }

  const Point& x1 = nodes_[0]->pt();
  const Point& x2 = nodes_[1]->pt();

  int nextGID = mesh->nextGID();
  midpoint_ = rcp(new Node(nextGID, 0.5*(x1 + x2), ownerProc_));
  mesh->addNode(midpoint_);

  int s;
  RCP<Edge> sub1 = mesh->tryEdge(nodes_[0], midpoint_, s);
  RCP<Edge> sub2 = mesh->tryEdge(midpoint_, nodes_[1], s);
  sub1->setParent(this);
  sub2->setParent(this);

  sub1->setLabel(label_);
  sub2->setLabel(label_);

  setChildren(sub1.get(), sub2.get());

  return midpoint_;
}

