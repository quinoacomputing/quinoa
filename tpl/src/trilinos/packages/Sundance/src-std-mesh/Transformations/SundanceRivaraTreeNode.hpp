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


#ifndef SUNDANCERIVARATREENODE_H
#define SUNDANCERIVARATREENODE_H

#include "SundanceDefs.hpp"

namespace Sundance
{
  namespace Rivara
    {
      /**
       * Class TreeNode represents a node in a Rivara mesh refinement tree.
       * Each node has either zero or two children depending on whether
       * it has been refined. All non-root nodes have a pointer back to their
       * parent.
       *
       * Only maximal elements will be responsible for deleting their children.
       * Therefore, the TreeNode dtor does not delete children; subtypes
       * that need to delete children should call the deleteChildren()
       * method.
       */
      class TreeNode 
        {
        public:
          /** Empty ctor */
          TreeNode();

          virtual ~TreeNode(){;}

          /** Delete the node's children */
          void deleteChildren();


          /** set the parent of this node */
          void setParent(TreeNode* parent) {parent_ = parent;}

          /** Set the two children of this node */
          void setChildren(TreeNode* left, TreeNode* right)
            {left_ = left; right_ = right;}

          /** return the leftmost leaf beneath this node */
          const TreeNode* first() const ;

          /** return the rightmost leaf beneath this node */
          const TreeNode* last() const ;

          /** Indicate whether this is the leftward child of another node */
          bool isLeftChild() const ;

          /** Indicate whether this is the rightward child of another node */
          bool isRightChild() const ;

          /** Return the next leaf in a left-to-right walk of the tree. If this
           * is the last leaf, return 0. */
          const TreeNode* next() const ;

          /** Indicate whether this node has children */
          bool hasChildren() const {return left_ != 0;}

          /** Return a count of the number of leaves */
          int numLeaves() const ;


        protected:
          TreeNode* left() {return left_;}

          TreeNode* right() {return right_;}
        private:

          TreeNode* parent_;

          TreeNode* left_;

          TreeNode* right_;
        };
    }
}

#endif
