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


#ifndef SUNDANCERIVARANODE_H
#define SUNDANCERIVARANODE_H

#include "SundanceDefs.hpp"
#include "SundancePoint.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"


namespace Sundance
{
namespace Rivara
{
class Element;
class Edge;
using Sundance::Point;
using Teuchos::Array;

/**
 * Class Node is a vertex in a simplicial mesh.
 */
class Node
{
public:
  Node(int gid, const Point& x, int ownerProc, int label=-1);

  /**
   * Return the spatial position of this node
   */
  const Point& pt() const ;

  /**
   * Return the local index of this node
   */
  int localIndex() const {return localIndex_;}


  /**
   * Return the global index of this node
   */
  int globalIndex() const {return globalIndex_;}

  /**
   * Set the local index of this node
   */
  void setLocalIndex(int localIndex) {localIndex_ = localIndex;}

  /**
   * Add an element to the list of elements containing this node.
   */
  void addConnectingElement(Element* elem);

  /**
   * Add an edge to the list of edges containing this node
   */
  void addConnectingEdge(Edge* edge);

  /**
   * Return the rank of the proc that owns this node
   */
  int ownerProc() const {return ownerProc_;}


  /**
   * Return the label of this node
   */
  int label() const {return label_;}


private:
  int label_;
  int localIndex_;
  int globalIndex_;
  Point x_;

  Array<Element*> elements_;
  Array<Edge*> edges_;

  int ownerProc_;
};
}
}



#endif
