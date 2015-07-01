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

#ifndef PLAYA_MULTI_VECTOR_OPERATOR_IMPL_HPP
#define PLAYA_MULTI_VECTOR_OPERATOR_IMPL_HPP

#include "PlayaMultiVectorOperatorDecl.hpp"
#include "PlayaVectorImpl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#endif

namespace Playa
{

/*
 * Construct from an array of vectors and a specifier for the 
 * domain space. 
 */
template <class Scalar> inline
MultiVectorOperator<Scalar>
::MultiVectorOperator(const Teuchos::Array<Vector<Scalar> >& cols,
  const VectorSpace<Scalar>& domain)
  : LinearOpWithSpaces<Scalar>(domain, cols[0].space()),
    cols_(cols)
{
  TEUCHOS_TEST_FOR_EXCEPTION(cols.size() == 0, std::runtime_error,
    "empty multivector given to MultiVectorOperator ctor");
  for (int i=1; i<cols.size(); i++)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(cols[i].space() != cols[0].space(), std::runtime_error,
      "inconsistent vector spaces in  MultiVectorOperator ctor");
  }
}


/*
 * Apply does an element-by-element multiply between the input 
 * vector, x, and the diagonal values.
 */
template <class Scalar> inline
void MultiVectorOperator<Scalar>
::apply(
  Teuchos::ETransp transApplyType,
  const Vector<Scalar>& in, 
  Vector<Scalar> out) const
{
  if (transApplyType == NO_TRANS)
  {
    out.zero();

    for (int i=0; i<cols_.size(); i++)
    {
      out.update(in[i], cols_[i]);
    }
  }
  else
  {
    out.zero();

    for (int i=0; i<cols_.size(); i++)
    {
      out[i] = in.dot(cols_[i]);
    }
  }
}



/* Return the kth row  */
template <class Scalar> inline
void MultiVectorOperator<Scalar>
::getRow(const int& k, 
  Teuchos::Array<int>& indices, 
  Teuchos::Array<Scalar>& values) const
{
  int low = this->range()->baseGlobalNaturalIndex();
  indices.resize(cols_.size());
  values.resize(cols_.size());
  for (int j=0; j<cols_.size(); j++)
  {
    indices[j] = j;
    values[j] = cols_[j][k-low];
  }
}


  
template <class Scalar> inline
LinearOperator<Scalar> multiVectorOperator(
  const Teuchos::Array<Vector<Scalar> >& cols,
  const VectorSpace<Scalar>& domain)
{
  RCP<LinearOperatorBase<Scalar> > A
    = rcp(new MultiVectorOperator<Scalar>(cols, domain));

  return A;
}


}

#endif
