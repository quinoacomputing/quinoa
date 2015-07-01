// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <assert.h>

#include <limits>

#include "ConstrainedOptPack_MatrixSymIdentitySerial.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace ConstrainedOptPack {

// Constructors

MatrixSymIdentitySerial::MatrixSymIdentitySerial(size_type size, value_type scale)
{
  this->initialize(size,scale);
}

void MatrixSymIdentitySerial::initialize(size_type size, value_type scale)
{
  size_  = size;
  scale_ = scale;
}

// Overridden from MatrixBase

size_type MatrixSymIdentitySerial::rows() const
{
  return size_;
}

size_type MatrixSymIdentitySerial::nz() const
{
  return size_;
}

// Overridden from MatrixOp

std::ostream& MatrixSymIdentitySerial::output(std::ostream& out) const
{
  out << "Identity matrix of size " << size_ << " x " << size_ << std::endl;
  return out;
}

// Overridden from MatrixOpSerial

void MatrixSymIdentitySerial::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
  ,const DVectorSlice& x, value_type b
  ) const
{
  DenseLinAlgPack::Vp_MtV_assert_sizes( y->dim(), rows(), cols(), BLAS_Cpp::no_trans, x.dim() );
  DenseLinAlgPack::Vt_S(y,b);
  DenseLinAlgPack::Vp_StV(y,a*scale_,x);
}

// Overridden from MatrixNonsinguarSerial

void MatrixSymIdentitySerial::V_InvMtV(
  DVectorSlice* y, BLAS_Cpp::Transp M_trans, const DVectorSlice& x
  ) const
{
  DenseLinAlgPack::Vp_MtV_assert_sizes( y->dim(), rows(), cols(), BLAS_Cpp::no_trans, x.dim() );
  LinAlgOpPack::V_StV(y,scale_,x);
}

// Overridden from MatrixSymNonsing

void MatrixSymIdentitySerial::M_StMtInvMtM(
    DMatrixSliceSym* S, value_type a
    ,const MatrixOpSerial& B, BLAS_Cpp::Transp B_trans
    ,EMatrixDummyArg dummy_arg
  ) const
{
  this->MatrixSymNonsingSerial::M_StMtInvMtM(S,a,B,B_trans,dummy_arg);
  // ToDo: Implement by calling S = b*S + scale*a*op(B')*op(B)
}

// Overridden from MatrixExtractInvCholFactor

void MatrixSymIdentitySerial::extract_inv_chol( DMatrixSliceTriEle* InvChol ) const
{
  if( scale_ < 0.0 )
    throw std::logic_error(
      "MatrixSymIdentitySerial::extract_inv_chol(...) : "
      "Error, we can not compute the inverse cholesky factor "
      "of a negative definite matrix." );
  DenseLinAlgPack::assign( &InvChol->gms(), 0.0 );
  InvChol->gms().diag() = 1.0 / std::sqrt(scale_);
}

} // end namespace ConstrainedOptPack
