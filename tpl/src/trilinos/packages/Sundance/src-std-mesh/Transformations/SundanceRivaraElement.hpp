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


#ifndef SUNDANCERIVARAELEMENT_H
#define SUNDANCERIVARAELEMENT_H

#include "SundanceDefs.hpp"
#include "SundanceRivaraTreeNode.hpp"
#include "SundanceRivaraNode.hpp"
#include "SundanceRivaraEdge.hpp"
#include "SundanceRivaraFace.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace Sundance
{
  namespace Rivara
  {
    class RivaraMesh;
    using Teuchos::Array;
    using Teuchos::RefCountPtr;

    /**
     * class Element is a maximal-dimension element in a simplicial mesh.
     */
    class Element : public TreeNode
    {
    public:
      Element(RivaraMesh* mesh,
        const RCP<Node>& a,
        const RCP<Node>& b,
        const RCP<Node>& c,
        int ownerProc, int label);

      Element(RivaraMesh* mesh,
        const RCP<Node>& a,
        const RCP<Node>& b,
        const RCP<Node>& c,
        const RCP<Node>& d,
        int ownerProc, int label);


      /** dtor */
      virtual ~Element() {deleteChildren();}

      /** compute the volume of this element. */
      double volume() const ;

      /** Indicate whether any nodes are hanging at this point
       * in the refinement algorithm.
       */
      bool hasHangingNode() const ;

      /**
       * Return the rank of the proc that owns this element.
       */
      int ownerProc() const {return ownerProc_;}

      /**
       * Return the label
       */
      int label() const {return label_;}

      /**
       * Return the index of the longest edge
       */
      int longestEdgeIndex() const ;

      /**
       * Refine the element
       */
      void refine(RivaraMesh* mesh, double maxArea);

      /**
       * Return a list of the neighbors of this element
       */
      void getNeighbors(Array<Element*>& neighbors,
        Array<int>& weights) const ;

      /**
       * Return the element's nodes
       */
      const Array<RCP<Node> >& nodes() const
      {return nodes_;}

      /** */
      const RCP<Edge>& edge(int i) const {return edges_[i];}
      /** */
      RCP<Edge> edge(int i) {return edges_[i];}
      
      /** */
      const RCP<Face>& face(int i) const {return faces_[i];}
      /** */
      RCP<Face> face(int i) {return faces_[i];}

      /** */
      bool hasNoEdgeLabels() const ;

      /** */
      bool hasNoFaceLabels() const ;

      /**
       * Return the element's nodes
       */
      Array<int> showNodes() const ;
    private:

      int label_;

      Array<RCP<Node> > nodes_;

      Array<RCP<Edge> > edges_;

      Array<RCP<Face> > faces_;

      Array<int> edgeSigns_;

      int ownerProc_;
    };

  }
}

#endif
