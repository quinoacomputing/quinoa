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

#ifndef PLAYA_SIMPLE_ADDED_OP_IMPL_HPP
#define PLAYA_SIMPLE_ADDED_OP_IMPL_HPP



#include "PlayaSimpleAddedOpDecl.hpp"
#include "PlayaSimpleZeroOpDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_Array.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaSimpleZeroOpImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;



/*
 * Represent a sum of operators A_0 + A_1 + ... + A_n.
 */
template <class Scalar> inline
SimpleAddedOp<Scalar>::SimpleAddedOp(
  const Array<LinearOperator<Scalar> >& ops)
  : LinearOpWithSpaces<Scalar>(
    ops[0].domain(), ops[0].range()
    ) 
  , ops_(ops)
{
  TEUCHOS_TEST_FOR_EXCEPT(ops_.size() <= 1);
  for (int i=1; i<ops_.size(); i++)
  {
    TEUCHOS_TEST_FOR_EXCEPT(!(ops[i].range() == ops[0].range()));
    TEUCHOS_TEST_FOR_EXCEPT(!(ops[i].domain() == ops[0].domain()));
  }
}
  
/* */
template <class Scalar> inline
void SimpleAddedOp<Scalar>::apply(Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in,
  Vector<Scalar> out) const
{
  Tabs tab(0);
  PLAYA_MSG2(this->verb(), tab << "SimpleAddedOp::apply()");

  Vector<Scalar> tmp=out.copy();
  tmp.zero();
  for (int i=0; i<ops_.size(); i++)
  {
    Tabs tab1;
    PLAYA_MSG3(this->verb(), tab1 << "applying term i=" << i << " of " 
      << ops_.size());
    if (transApplyType == Teuchos::NO_TRANS)
      tmp += ops_[i] * in;
    else if (transApplyType == Teuchos::TRANS)
      tmp += ops_[i].transpose() * in;
    else 
      TEUCHOS_TEST_FOR_EXCEPT(transApplyType != Teuchos::TRANS && transApplyType != Teuchos::NO_TRANS);
  }
  out.acceptCopyOf(tmp);

  PLAYA_MSG2(this->verb(), tab << "done SimpleAddedOp::apply()");
}
  
/* */
template <class Scalar> inline
std::string SimpleAddedOp<Scalar>::description() const 
{
  std::string rtn="(";
  for (int i=0; i<ops_.size(); i++)
  {
    if (i > 0) rtn += "+";
    rtn += ops_[i].description();
  }
  rtn += ")";
  return rtn;
}



template <class Scalar> inline
LinearOperator<Scalar> addedOperator(
  const Array<LinearOperator<Scalar> >& ops)
{
  /* We will strip out any zero operators */
  Array<LinearOperator<Scalar> > strippedOps;

  for (int i=0; i<ops.size(); i++)
  {
    LinearOperator<Scalar> op_i = ops[i];

    /* Ignore any zero operators */
    const SimpleZeroOp<Scalar>* zPtr 
      = dynamic_cast<const SimpleZeroOp<Scalar>*>(op_i.ptr().get());

    if (zPtr != 0) continue;

    strippedOps.append(op_i);
  }
  
  TEUCHOS_TEST_FOR_EXCEPT(strippedOps.size() < 1);
  if (strippedOps.size()==1) return strippedOps[0];
  
  RCP<LinearOperatorBase<Scalar> > op 
    = rcp(new SimpleAddedOp<Scalar>(strippedOps));
  
  return op;
}

template <class Scalar> inline
LinearOperator<Scalar> operator+(const LinearOperator<Scalar>& A,
  const LinearOperator<Scalar>& B)
{
  return addedOperator(Array<LinearOperator<Scalar> >(tuple(A, B)));
}

}

#endif
