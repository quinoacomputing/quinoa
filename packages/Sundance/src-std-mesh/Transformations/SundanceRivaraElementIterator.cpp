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


#include "SundanceRivaraElementIterator.hpp"
#include "SundanceRivaraMesh.hpp"
#include "SundanceRivaraElement.hpp"
#include "SundanceRivaraNode.hpp"


using namespace Sundance::Rivara;
using namespace Teuchos;


ElementIterator::ElementIterator(const RivaraMesh* mesh)
	: mesh_(const_cast<RivaraMesh*>(mesh)), rootIndex_(-1), current_(0),
		startingNewTree_(true)
{}

bool ElementIterator::hasMoreElements() const 
{
	if (rootIndex_ < mesh_->numRootElements()-1)
		{
			return true;
		}
	else if (current_->next() != 0)
		{
			return true;
		}
	else return false;
}

const Element* ElementIterator::getNextElement() const 
{
	if (current_==0)
		{
			current_ = const_cast<Element*>(dynamic_cast<const Element*>(mesh_->element(0)->first()));
			rootIndex_ = 0;
		}
	else if (current_->next() != 0)
		{
			current_ = const_cast<Element*>(dynamic_cast<const Element*>(current_->next()));
		}
	else
		{
			rootIndex_++;
			current_ = const_cast<Element*>(dynamic_cast<const Element*>(mesh_->element(rootIndex_)->first()));
		}
	
	return current_;
	
}

	
