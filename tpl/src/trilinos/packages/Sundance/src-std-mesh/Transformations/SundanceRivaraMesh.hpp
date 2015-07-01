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


#ifndef SUNDANCERIVARAMESH_HPP
#define SUNDANCERIVARAMESH_HPP

#include "SundanceDefs.hpp"
#include "SundanceRivaraElement.hpp"
#include <stack>
#include "SundanceMap.hpp"
#include "SundanceIncrementallyCreatableMesh.hpp"
#include "SundanceRivaraElementIterator.hpp"


namespace Sundance
{
namespace Rivara
{
using Sundance::Map;


class RivaraMesh 
{
public:
  RivaraMesh(int dim, const MPIComm& comm);

  int addNode(const RCP<Node>& node);
  int addVertex(int globalIndex, const Point& x, int ownerProcID, int label);

  void addElement(const RCP<Element>& tri);
  int addElement(int globalIndex, const Array<int>& vertexGIDs, int ownerProc,
    int label);

  RCP<Edge> tryEdge(const RCP<Node>& a,
    const RCP<Node>& b,
    int& edgeSign);

  RCP<Face> tryFace(const RCP<Node>& a,
    const RCP<Node>& b,
    const RCP<Node>& c);

  const RCP<Face>& getFace(const RCP<Node>& a,
    const RCP<Node>& b,
    const RCP<Node>& c) const ;

  const RCP<Node>& node(int i) const {return nodes_[i];}

  int numNodes() const {return nodes_.length();}

  std::stack<Element*>& refinementSet()
    {return refinementSet_;}

  std::stack<double>& refinementAreas()
    {return refinementAreas_;}

  void refine();

  ElementIterator iterator() const ;

  friend class ElementIterator;

  RCP<Element> element(int i) const {return elements_[i];}

  int numRootElements() const {return elements_.length();}

  int numElements() const ;

  int spatialDim() const ;

  int& nextGID() {return nextGID_;}

  int nextGID() const {return nextGID_;}
private:

  int spatialDim_;
  
  int nextGID_;

  Array<RCP<Node> > nodes_;

  Array<RCP<Edge> > edges_;

  Array<RCP<Face> > faces_;

  Array<RCP<Element> > elements_;

  Array<Map<int, int> > nodeToEdgeMap_;

  Map<FaceNodes, int> faceToLIDMap_;

  std::stack<Element*> refinementSet_;

  std::stack<double> refinementAreas_;
};
}
}

#endif
