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

#ifndef PLAYA_SIMPLE_SCALED_OP_IMPL_HPP
#define PLAYA_SIMPLE_SCALED_OP_IMPL_HPP



#include "PlayaSimpleScaledOpDecl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;





/*
 * --- scaled op
 */

template <class Scalar> inline
SimpleScaledOp<Scalar>::SimpleScaledOp(const Scalar& alpha,
  const LinearOperator<Scalar>& A)
  : LinearOpWithSpaces<Scalar>(
    A.domain(), A.range()
    ) 
  , alpha_(alpha), A_(A)
{}
  
/* */
template <class Scalar> inline
void SimpleScaledOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleScaledOp::apply()");

  if (transApplyType == Teuchos::NO_TRANS)
    A_.apply(in, out);
  else if (transApplyType == Teuchos::TRANS)
    A_.applyTranspose(in, out);
  else 
    TEUCHOS_TEST_FOR_EXCEPT(transApplyType !=Teuchos::TRANS && transApplyType != Teuchos::NO_TRANS);

  out.scale(alpha_);

  PLAYA_MSG2(this->verb(), tab << "done SimpleScaledOp::apply()");
}
  
/* */
template <class Scalar> inline
std::string SimpleScaledOp<Scalar>::description() const 
{
  return "ScaledOp[alpha="  + Teuchos::toString(alpha_)
    + ", " + A_.description() + "]";
}


/* */
template <class Scalar> inline
void SimpleScaledOp<Scalar>::print(std::ostream& os) const 
{
  Tabs tab(0);
  os << tab << "ScaledOp[" << std::endl;
  Tabs tab1;
  os << tab1 << "scale = " << alpha_ << std::endl;
  os << tab1 << "operator = " << A_.description() << std::endl;
  os << tab << "]" << std::endl;
}



template <class Scalar> inline
LinearOperator<Scalar> scaledOperator(
  const Scalar& scale,
  const LinearOperator<Scalar>& op)
{
  RCP<LinearOperatorBase<Scalar> > A 
    = rcp(new SimpleScaledOp<Scalar>(scale, op));

  return A;
}


template <class Scalar> inline
LinearOperator<Scalar> operator*(const Scalar& a, const LinearOperator<Scalar>& A)
{
  return scaledOperator(a, A);
}
  

}

#endif
