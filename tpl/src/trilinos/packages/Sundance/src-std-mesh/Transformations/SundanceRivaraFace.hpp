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


#ifndef SUNDANCERIVARAFACE_H
#define SUNDANCERIVARAFACE_H

#include "SundanceDefs.hpp"
#include "SundanceRivaraTreeNode.hpp"
#include "SundanceRivaraNode.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;
namespace Sundance
{
namespace Rivara
{
class Element;
class RivaraMesh;

using Sundance::Set;
using Sundance::makeSet;
using Teuchos::Array;
using Teuchos::RefCountPtr;

/** */
class FaceNodes
{
public:
  /** */
  FaceNodes(
    const RCP<Node>& a, 
    const RCP<Node>& b,
    const RCP<Node>& c)
    : nodeLIDSet_(rcp(new Set<int>(makeSet(a->localIndex(), b->localIndex(), c->localIndex()))))
    {}


  /** */
  bool operator<(const FaceNodes& other) const
    {
      return (*nodeLIDSet_) < (*(other.nodeLIDSet_));
    }

  /** */
  const Set<int>& nodes() const {return *nodeLIDSet_;}

private:
  RCP<Set<int> > nodeLIDSet_;
};

/**
 * class Face is a two-dimensional face in a simplicial mesh.
 */

class Face
{
public:
  /**
   * Construct with three nodes
   */
  Face(
    const RCP<Node>& a, 
    const RCP<Node>& b,
    const RCP<Node>& c
    )
    : label_(-1), 
      globalIndex_(-1), 
      nodes_(a,b,c), 
      ownerProc_(std::max(std::max(a->ownerProc(), b->ownerProc()), c->ownerProc()))
    {}

  /** */
  const FaceNodes& nodes() const {return nodes_;}

  /** */
  int ownerProc() const {return ownerProc_;}

  /**
   * Return the global index of this edge
   */
  int globalIndex() const {return globalIndex_;}

  /**
   * Set the global index of this edge
   */
  void setGlobalIndex(int globalIndex) {globalIndex_ = globalIndex;}

  /**
   * Set the label of this edge
   */
  void setLabel(int label) {label_=label;}

  /** 
   * Get the label
   */
  int label() const {return label_;}
private:

  int label_;
  int globalIndex_;
  FaceNodes nodes_;
  int ownerProc_;
};
}
}

#endif
