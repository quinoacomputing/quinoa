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

#ifndef PLAYA_SIMPLEBLOCKOP_IMPL_HPP
#define PLAYA_SIMPLEBLOCKOP_IMPL_HPP

#include "PlayaTabs.hpp"
#include "PlayaExceptions.hpp"
#include "PlayaSimpleBlockOpDecl.hpp"
#include "PlayaSimpleZeroOpDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleZeroOpImpl.hpp"
#endif



namespace Playa
{
using namespace Teuchos;



/* ---- Simplified linear op with spaces ------- */

template <class Scalar> inline
SimpleBlockOp<Scalar>::SimpleBlockOp(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range)
  : LinearOpWithSpaces<Scalar>(domain, range), blocks_(range.numBlocks())
{
  for (int i=0; i<blocks_.size(); i++) 
  {
    blocks_[i] = Array<LinearOperator<Scalar> >(domain.numBlocks());
    for (int j=0; j<blocks_[i].size(); j++)
    {
      blocks_[i][j] = zeroOperator(domain.getBlock(j), range.getBlock(i));
    }
  }
}

template <class Scalar> inline
int SimpleBlockOp<Scalar>::numBlockRows() const
{
  return blocks_.size();
}

template <class Scalar> inline
int SimpleBlockOp<Scalar>::numBlockCols() const
{
  return blocks_[0].size();
}

template <class Scalar> inline
const LinearOperator<Scalar>& SimpleBlockOp<Scalar>::getBlock(int i, int j) const 
{
  return blocks_[i][j];
}

template <class Scalar> inline
LinearOperator<Scalar> SimpleBlockOp<Scalar>::getNonconstBlock(int i, int j) 
{
  return blocks_[i][j];
}


template <class Scalar> inline
void SimpleBlockOp<Scalar>::setBlock(int i, int j, 
  const LinearOperator<Scalar>& Aij) 
{
  blocks_[i][j] = Aij;
}

template <class Scalar> inline
void SimpleBlockOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleBlockOp::apply()");
  if (transApplyType == Teuchos::NO_TRANS)
  {
    out.zero();
    for (int i=0; i<this->numBlockRows(); i++)
    {
      for (int j=0; j<this->numBlockCols(); j++)
      {
        Vector<Scalar> tmp; 
        blocks_[i][j].apply(in.getBlock(j), tmp);
        out.getNonConstBlock(i).update(1.0, tmp);
      }
    }
  }

  else if (transApplyType == Teuchos::TRANS)
  {
    for (int i=0; i<this->numBlockCols(); i++)
    {
      out.zero();
      for (int j=0; j<this->numBlockRows(); j++)
      {
        Vector<Scalar> tmp;
        blocks_[j][i].applyTranspose(in.getBlock(j),tmp);
        out.getNonConstBlock(i).update(1.0, tmp);
      }
    }
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPT(transApplyType != Teuchos::TRANS && transApplyType != Teuchos::NO_TRANS);
  }
  PLAYA_MSG2(this->verb(), tab << "done SimpleBlockOp::apply()");
}


template <class Scalar> inline
LinearOperator<Scalar> makeBlockOperator(
  const VectorSpace<Scalar>& domain,
  const VectorSpace<Scalar>& range
  )
{
  RCP<SimpleBlockOp<Scalar> > b = 
    rcp(new SimpleBlockOp<Scalar>(domain, range));
  RCP<LinearOperatorBase<Scalar> > p = b;
  return p;
}



}

#endif
