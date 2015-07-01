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

#include "ConstrainedOptPack_MatrixIdentConcatStd.hpp"
#include "Teuchos_Assert.hpp"

namespace ConstrainedOptPack {

// Setup and representation access

MatrixIdentConcatStd::MatrixIdentConcatStd()
{
  this->set_uninitialized();
}

void MatrixIdentConcatStd::initialize(
    const VectorSpace::space_ptr_t&    space_cols
    ,const VectorSpace::space_ptr_t&   space_rows
    ,ETopBottom                        top_or_bottom
    ,value_type                        alpha
    ,const D_ptr_t                     &D_ptr
    ,BLAS_Cpp::Transp                  D_trans
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    space_cols.get() == NULL, std::invalid_argument
    ,"MatrixIdentConcatStd::initialize(...): Error, "
    "space_cols.get() can not be NULL!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    space_rows.get() == NULL, std::invalid_argument
    ,"MatrixIdentConcatStd::initialize(...): Error, "
    "space_rows.get() can not be NULL!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    D_ptr.get() == NULL, std::invalid_argument
    ,"MatrixIdentConcatStd::initialize(...): Error, "
    "D_ptr.get() can not be NULL!" );
#endif
  const size_type
    D_rows   = D_ptr->rows(),
    D_cols   = D_ptr->cols(),
    opD_rows = BLAS_Cpp::rows( D_rows, D_cols, D_trans ),
    opD_cols = BLAS_Cpp::cols( D_rows, D_cols, D_trans ),
    rows     = opD_rows + opD_cols;
  space_cols_ = space_cols;
  space_rows_ = space_rows;
  alpha_      = alpha;
  D_ptr_      = D_ptr;
  D_trans_    = D_trans;
  D_rng_      = top_or_bottom == TOP ? Range1D(1,opD_rows)      : Range1D(opD_cols+1,rows);
  I_rng_      = top_or_bottom == TOP ? Range1D(opD_rows+1,rows) : Range1D(1,opD_cols);
}

void MatrixIdentConcatStd::set_uninitialized()
{
  namespace rcp = MemMngPack;
  space_cols_ = Teuchos::null;
  space_rows_ = Teuchos::null;
  alpha_      = 0.0;
  D_ptr_      = Teuchos::null;
  D_trans_    = BLAS_Cpp::no_trans;
  D_rng_      = Range1D::Invalid;
  I_rng_      = Range1D::Invalid;
}

const MatrixIdentConcatStd::D_ptr_t& MatrixIdentConcatStd::D_ptr() const
{
  return D_ptr_;
}

// Overridden form MatrixIdentConcat

Range1D MatrixIdentConcatStd::D_rng() const
{
  return D_rng_;
}

Range1D MatrixIdentConcatStd::I_rng() const
{
  return I_rng_;
}

value_type MatrixIdentConcatStd::alpha() const
{
  return alpha_;
}

const MatrixOp& MatrixIdentConcatStd::D() const
{
  return *D_ptr_;
}

BLAS_Cpp::Transp MatrixIdentConcatStd::D_trans() const
{
  return D_trans_;
}

// Overridden from MatrixOp

const VectorSpace& MatrixIdentConcatStd::space_cols() const
{
  return *space_cols_;
}

const VectorSpace& MatrixIdentConcatStd::space_rows() const
{
  return *space_rows_;
}

MatrixOp& MatrixIdentConcatStd::operator=(const MatrixOp& m)
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // Finish!
  return *this;
}

// private

void MatrixIdentConcatStd::assert_initialized() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    space_cols_.get() == NULL, std::logic_error
    ,"Error, the MatrixIdentConcatStd object has not been initialized!" );
}

} // end namespace ConstrainedOptPack
