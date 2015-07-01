/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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

#ifndef PLAYA_BLOCKEDVECTORSPACEBASEIMPL_HPP
#define PLAYA_BLOCKEDVECTORSPACEBASEIMPL_HPP

#include "PlayaBlockVectorSpaceDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"

namespace Playa
{
template <class Scalar> inline
int BlockVectorSpaceBase<Scalar>::dim() const 
{
  int rtn = 0;
  for (int i=0; i<this->numBlocks(); i++) 
  {
    rtn += this->getBlock(i).dim();
  }
  return rtn;
}

template <class Scalar> inline
int BlockVectorSpaceBase<Scalar>::numLocalElements() const 
{
  int rtn = 0;
  for (int i=0; i<this->numBlocks(); i++) 
  {
    rtn += this->getBlock(i).numLocalElements();
  }
  return rtn;
}

template <class Scalar> inline
bool BlockVectorSpaceBase<Scalar>::isCompatible(const VectorSpaceBase<Scalar>* other) const 
{
  if (this->numBlocks() != other->numBlocks()) return false;

  const BlockVectorSpaceBase<Scalar>* bs 
    = dynamic_cast<const BlockVectorSpaceBase<Scalar>* >(other);
    
  if (bs == 0) return false;

  for (int i=0; i<this->numBlocks(); i++) 
  {
    if (! (this->getBlock(i).isCompatible(bs->getBlock(i))))
      return false;
  }
  return true;
}


template <class Scalar> inline
std::string BlockVectorSpaceBase<Scalar>::description() const 
{
  std::ostringstream rtn;
  rtn << "BlockVS[";
  for (int b=0; b<this->numBlocks(); b++)
  {
    if (b > 0) rtn << ", ";
    rtn << this->getBlock(b).description();
  }
  rtn << "]";
  return rtn.str();
}


  
  
}

#endif
