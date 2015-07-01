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

#include <ostream>

#include "AbstractLinAlgPack_MatrixSymIdent.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

namespace AbstractLinAlgPack {

// Constructors/initalizers

MatrixSymIdent::MatrixSymIdent(
  const VectorSpace::space_ptr_t&          vec_space
  ,const value_type                        scale
  )
{
  this->initialize(vec_space,scale);
}

void MatrixSymIdent::initialize(
  const VectorSpace::space_ptr_t&          vec_space
  ,const value_type                        scale
  )
{
  vec_space_ = vec_space;
  scale_     = scale;
}

// Overridden from MatrixBase

size_type MatrixSymIdent::rows() const
{
  return vec_space_.get() ? vec_space_->dim() : 0;
}

size_type MatrixSymIdent::nz() const
{
  return vec_space_.get() ? vec_space_->dim() : 0;
}

// Overridden from MatrixOp

const VectorSpace& MatrixSymIdent::space_cols() const {
  return *vec_space_;
}

std::ostream& MatrixSymIdent::output(std::ostream& out) const
{
  out << "Identity matrix of dimension " << rows() << " x " << rows() << std::endl;
  return out;
}

void MatrixSymIdent::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  ,const Vector& x, value_type b
  ) const
{
  AbstractLinAlgPack::Vp_MtV_assert_compatibility( y, *this, BLAS_Cpp::no_trans, x );
  Vt_S(y,b);
    Vp_StV(y,a*scale_,x);
}

// Overridden from MatrixNonsing

void MatrixSymIdent::V_InvMtV(
  VectorMutable* y, BLAS_Cpp::Transp M_trans, const Vector& x
  ) const
{
  AbstractLinAlgPack::Vp_MtV_assert_compatibility( y, *this, BLAS_Cpp::no_trans, x );
  LinAlgOpPack::V_StV(y,1.0/scale_,x);
}

} // end namespace AbstractLinAlgPack
