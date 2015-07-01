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

#ifndef Playa_MULTI_VECTOR_OPERATOR_DECL_HPP
#define Playa_MULTI_VECTOR_OPERATOR_DECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaRowAccessibleOp.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaHandleable.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "PlayaVectorDecl.hpp"

namespace Playa
{
/** 
 * A MultiVectorOperator is a linear operator whose
 * rows or columns are represented as a multivector 
 */
template <class Scalar> 
class MultiVectorOperator 
  : public LinearOpWithSpaces<Scalar>,
    public RowAccessibleOp<Scalar>
{
public:

  /**
   * Construct from an array of vectors and a specifier for the 
   * domain space. 
   */
  MultiVectorOperator(const Teuchos::Array<Vector<Scalar> >& cols,
    const VectorSpace<Scalar>& domain);

  /** Virtual dtor */
  virtual ~MultiVectorOperator(){;}


  /** 
   * Apply does an element-by-element multiply between the input 
   * vector, x, and the diagonal values.
   */
  virtual void apply(
    Teuchos::ETransp transType,
    const Vector<Scalar>& in, 
    Vector<Scalar> out) const ;


  /** Return the kth row  */
  void getRow(const int& k, 
    Teuchos::Array<int>& indices, 
    Teuchos::Array<Scalar>& values) const ;

private:

  Teuchos::Array<Vector<Scalar> > cols_;
};

/** \relates MultiVectorOperator */
template <class Scalar> 
LinearOperator<Scalar> multiVectorOperator(
  const Teuchos::Array<Vector<Scalar> >& cols,
  const VectorSpace<Scalar>& domain);


}

#endif
