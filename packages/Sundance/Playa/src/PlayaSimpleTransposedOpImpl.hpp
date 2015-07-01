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

#ifndef PLAYA_SIMPLE_TRANSPOSED_OP_IMPL_HPP
#define PLAYA_SIMPLE_TRANSPOSED_OP_IMPL_HPP



#include "PlayaSimpleTransposedOpDecl.hpp"
#include "PlayaSimpleZeroOpDecl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleZeroOpImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;




/*
 * --- transposed op
 */

template <class Scalar> inline
SimpleTransposedOp<Scalar>::SimpleTransposedOp(const LinearOperator<Scalar>& A)
  : LinearOpWithSpaces<Scalar>(
    A.range(), A.domain()
    ) 
  , A_(A)
{}
  
/* */
template <class Scalar> inline
void SimpleTransposedOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleTransposedOp::apply()");

  if (transApplyType == Teuchos::NO_TRANS)
    A_.applyTranspose(in, out);
  else if (transApplyType == Teuchos::TRANS)
    A_.apply(in, out);
  else 
    TEUCHOS_TEST_FOR_EXCEPT(transApplyType !=Teuchos::TRANS && transApplyType != Teuchos::NO_TRANS);

  PLAYA_MSG2(this->verb(), tab << "done SimpleTransposedOp::apply()");
}
  
/* */
template <class Scalar> inline
std::string SimpleTransposedOp<Scalar>::description() const 
{
  return "(" + A_.description() + "^T)";
}



template <class Scalar> inline
LinearOperator<Scalar> transposedOperator(
  const LinearOperator<Scalar>& op)
{

  /* If the operator is a transpose, return the untransposed op */
  const SimpleTransposedOp<Scalar>* tPtr
    = dynamic_cast<const SimpleTransposedOp<Scalar>*>(op.ptr().get());
  if (tPtr)
  {
    return tPtr->op();
  }

  /* If the operator is zero, return a transposed zero */
  const SimpleZeroOp<Scalar>* zPtr 
    = dynamic_cast<const SimpleZeroOp<Scalar>*>(op.ptr().get());

  if (zPtr != 0) 
  {
    VectorSpace<Scalar> r = op.range();
    VectorSpace<Scalar> d = op.domain();
    return zeroOperator(r, d);
  }


  /* Return a transposed operator */
  RCP<LinearOperatorBase<Scalar> > A
    = rcp(new SimpleTransposedOp<Scalar>(op));
      
  return A;
}

  


}

#endif
