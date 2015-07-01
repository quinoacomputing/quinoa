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


#include "SundanceRivaraTreeNode.hpp"
#include "PlayaExceptions.hpp"

using namespace Sundance::Rivara;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

TreeNode::TreeNode()
  : parent_(0), left_(0), right_(0)
{;}

void TreeNode::deleteChildren()
{
  if (left_ != 0) delete left_;
  if (right_ != 0) delete right_;
}

const TreeNode* TreeNode::first() const
{
  /* keep walking left as long as possible */
  if (left_ != 0) return left_->first();
  return this;
}

const TreeNode* TreeNode::last() const
{
  /* keep walking right as long as possible */
  if (right_ != 0) return right_->last();
  return this;
}

bool TreeNode::isRightChild() const
{
  /* return true if I am the right child of my parent */
  if (parent_!=0 && parent_->right_ == this) return true;
  return false;
}

bool TreeNode::isLeftChild() const
{
  /* return true if I am the left child of my parent */
  if (parent_!=0 && parent_->left_ == this) return true;
  return false;
}

const TreeNode* TreeNode::next() const 
{

  /* walk up the tree until we are at either the root or are the
   * left child */
  TreeNode* pos = const_cast<TreeNode*>(this);

  while (pos->isRightChild())
    {
      pos = pos->parent_;
    }

  /* if we are the left child, begin at the leftmost leaf of
   * the right sibling tree */
  if (pos->isLeftChild())
    {
      return pos->parent_->right_->first();
    }

  /* if we are at the root, there are no more unwalked leaves */
  if (pos->parent_==0) return pos->parent_;

  /* if we get to this point, there's a bug in this code */
  TEUCHOS_TEST_FOR_EXCEPT(true);

  return pos->parent_;
}

int TreeNode::numLeaves() const
{
  if (hasChildren()) 
		{
			return left_->numLeaves() + right_->numLeaves();
		}
  return 1;
}





