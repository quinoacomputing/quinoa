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


#ifndef SUNDANCERIVARAEDGE_H
#define SUNDANCERIVARAEDGE_H

#include "SundanceDefs.hpp"
#include "SundanceRivaraTreeNode.hpp"
#include "SundanceRivaraNode.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Sundance
{
  namespace Rivara
  {
    class Element;
    class RivaraMesh;

    using Teuchos::Array;
    using Teuchos::RefCountPtr;

    /**
     * class Edge is a one-dimensional edge in a simplicial mesh.
     */

    class Edge : public TreeNode
    {
    public:
      /**
       * Construct with two nodes
       */
      Edge(const RCP<Node>& a,
        const RCP<Node>& b);

      /**
       * Add an element to the list of elements containing this edge
       */
      void addConnectingElement(Element* elem);

      /**
       * Return the length of the edge.
       */
      double length() const ;

      /**
       * Return a list of the cofacets of this edge that still need refinement
       */
      void getUnrefinedCofacets(Array<Element*>& elements) const ;

      /**
       * Bisect the edge.
       * @return a new node created at the midpoint of the edge
       */
      RCP<Node> bisect(RivaraMesh* mesh);

      const RCP<Node>& node(int i) const {return nodes_[i];}

      int ownerProc() const {return ownerProc_;}

      /**
       * Return the global index of this edge
       */
      int globalIndex() const ;

      /**
       * Set the global index of this edge
       */
      void setGlobalIndex(int globalIndex) ;

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
      Array<RCP<Node> > nodes_;
      Array<Element*> elements_;

      RCP<Node> midpoint_;

      int ownerProc_;
    };
  }
}

#endif
