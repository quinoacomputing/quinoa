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

#ifndef PLAYA_DEFAULTBLOCKVECTORSPACEIMPL_HPP
#define PLAYA_DEFAULTBLOCKVECTORSPACEIMPL_HPP

#include "PlayaDefaultBlockVectorSpaceDecl.hpp"
#include "PlayaDefaultBlockVectorDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaExceptions.hpp"

namespace Playa
{

template <class Scalar> inline
DefaultBlockVectorSpace<Scalar>::
DefaultBlockVectorSpace(const Array<VectorSpace<Scalar> >& blocks) 
  : blocks_(blocks), baseGNI_(-1) 
{
  baseGNI_ = this->accumulateBaseGNI();
}

template <class Scalar> inline
const VectorSpace<Scalar>&  
DefaultBlockVectorSpace<Scalar>::getBlock(int b) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(b < 0 || b >= this->numBlocks(),
    RuntimeError, "block index b=" << b << " into vector space "
    << this->description() << " out of range [0,"
    << this->numBlocks() << ")");
  return blocks_[b];
}


template <class Scalar> inline
RCP<VectorBase<Scalar> > DefaultBlockVectorSpace<Scalar>
::createMember(const VectorSpace<Scalar>& self) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(this != self.ptr().get(), RuntimeError,
    "inconsistent self-reference in DefaultBlockVectorSpace::"
    "createMember()");
  return rcp(new DefaultBlockVector<Scalar>(self));
}


template <class Scalar> inline
const MPIComm& DefaultBlockVectorSpace<Scalar>::comm() const
{
  return getBlock(0).comm();
}

template <class Scalar> inline
int DefaultBlockVectorSpace<Scalar>::baseGlobalNaturalIndex() const
{
  return baseGNI_;
}




/** \relates VectorSpace */
template <class Scalar> inline
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1)
{
  Array<VectorSpace<Scalar> > x(1);
  x[0] = v1;
  return blockSpace<Scalar>(x);
}


/** \relates VectorSpace */
template <class Scalar> inline
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1,
  const VectorSpace<Scalar>& v2)
{
  Array<VectorSpace<Scalar> > x(2);
  x[0] = v1;
  x[1] = v2;
  return blockSpace<Scalar>(x);
}

/** \relates VectorSpace */
template <class Scalar> inline
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1,
  const VectorSpace<Scalar>& v2,
  const VectorSpace<Scalar>& v3)
{
  Array<VectorSpace<Scalar> > x(3);
  x[0] = v1;
  x[1] = v2;
  x[2] = v3;
  return blockSpace<Scalar>(x);
}

/** \relates VectorSpace */
template <class Scalar> inline
VectorSpace<Scalar> blockSpace(
  const VectorSpace<Scalar>& v1,
  const VectorSpace<Scalar>& v2,
  const VectorSpace<Scalar>& v3,
  const VectorSpace<Scalar>& v4)
{
  Array<VectorSpace<Scalar> > x(4);
  x[0] = v1;
  x[1] = v2;
  x[2] = v3;
  x[3] = v4;
  return blockSpace<Scalar>(x);
}

/** \relates VectorSpace */
template <class Scalar> inline
VectorSpace<Scalar> blockSpace(const Array<VectorSpace<Scalar> >& x)
{
  RCP<const VectorSpaceBase<Scalar> > rtn 
    = rcp(new DefaultBlockVectorSpace<Scalar>(x));
  return rtn;
}


}


 

#endif
