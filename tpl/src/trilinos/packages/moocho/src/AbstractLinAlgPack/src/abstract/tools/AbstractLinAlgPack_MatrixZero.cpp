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

#include "AbstractLinAlgPack_MatrixZero.hpp"
#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "Teuchos_Assert.hpp"

namespace AbstractLinAlgPack {

// Constructors/initializers

MatrixZero::MatrixZero(
  const VectorSpace::space_ptr_t&    space_cols
  ,const VectorSpace::space_ptr_t&   space_rows
  )
{
  this->initialize(space_cols,space_rows);
}

void MatrixZero::initialize(
  const VectorSpace::space_ptr_t&    space_cols
  ,const VectorSpace::space_ptr_t&   space_rows
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    (space_cols.get() == NULL && space_rows.get() != NULL)
    || (space_cols.get() != NULL && space_rows.get() == NULL)
    , std::invalid_argument
    ,"MatrixZero::initialize(...) : Error, the space_cols.get() and "
    "space_rows.get() must both be != NULL or == NULL" );
  space_cols_ = space_cols;
  space_rows_ = space_rows;
}

// Overridden from MatrixBase

size_type MatrixZero::rows() const
{
  return space_cols_.get() ? space_cols_->dim() : 0;
}

size_type MatrixZero::cols() const
{
  return space_rows_.get() ? space_rows_->dim() : 0;
}

size_type MatrixZero::nz() const
{
  return 0;
}

// Overridden form MatrixOp

const VectorSpace& MatrixZero::space_cols() const
{
  assert_initialized();
  return *space_cols_;
}

const VectorSpace& MatrixZero::space_rows() const
{
  assert_initialized();
  return *space_rows_;
}

void MatrixZero::zero_out()
{
  assert_initialized();
  // Automatically satisfied!
}

void MatrixZero::Mt_S( value_type alpha )
{
  assert_initialized();
  // Automatically satisfied!
}

MatrixOp& MatrixZero::operator=(const MatrixOp& M)
{
  assert_initialized();
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
  return *this;
}

std::ostream& MatrixZero::output(std::ostream& out) const
{
  assert_initialized();
  return out << "Zero matrix of dimension " << rows() << " x " << cols() << std::endl;
}

// Level-1 BLAS

bool MatrixZero::Mp_StM(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs) const
{
  assert_initialized();
  return true; // Nothing to do!
}

bool MatrixZero::Mp_StMtP(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  ) const
{
  assert_initialized();
  return true; // Nothing to do!
}

bool MatrixZero::Mp_StPtM(
  MatrixOp* m_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  , BLAS_Cpp::Transp M_trans
  ) const
{
  assert_initialized();
  return true; // Nothing to do!
}

bool MatrixZero::Mp_StPtMtP(
  MatrixOp* m_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  ) const
{
  assert_initialized();
  return true; // Nothing to do!
}

// Level-2 BLAS

void MatrixZero::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans_in
  , const Vector& x, value_type b
  ) const
{
  assert_initialized();
  Vt_S(y,b);
}

void MatrixZero::Vp_StMtV(
  VectorMutable* y, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& x, value_type b) const
{
  assert_initialized();
  Vt_S(y,b);
}

void MatrixZero::Vp_StPtMtV(
  VectorMutable* y, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const Vector& x, value_type b) const
{
  assert_initialized();
  Vt_S(y,b);
}

void MatrixZero::Vp_StPtMtV(
  VectorMutable* y, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const SpVectorSlice& x, value_type b) const
{
  assert_initialized();
  Vt_S(y,b);
}

value_type MatrixZero::transVtMtV(
  const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const Vector& v_rhs3) const
{
  assert_initialized();
  return 0.0; // Nothing to do!
}

value_type MatrixZero::transVtMtV(
  const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const SpVectorSlice& sv_rhs3) const
{
  assert_initialized();
  return 0.0; // Nothing to do!
}

void MatrixZero::syr2k(
  BLAS_Cpp::Transp M_trans, value_type alpha
  , const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
  , const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
  , value_type beta, MatrixSymOp* sym_lhs ) const
{
  assert_initialized();
  sym_lhs->Mt_S(beta);
}

// Level-3 BLAS

bool MatrixZero::Mp_StMtM(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOp& mwo_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
  assert_initialized();
  m_lhs->Mt_S(beta);
  return true;
}

bool MatrixZero::Mp_StMtM(
  MatrixOp* m_lhs, value_type alpha
  , const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  , BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
  assert_initialized();
  m_lhs->Mt_S(beta);
  return true;
}

bool MatrixZero::syrk(
  BLAS_Cpp::Transp M_trans, value_type alpha
  , value_type beta, MatrixSymOp* sym_lhs ) const
{
  assert_initialized();
  sym_lhs->Mt_S(beta);
  return true;
}

// private

void MatrixZero::assert_initialized() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    space_cols_.get() == NULL, std::logic_error
    ,"Error, the MatrixZero object has not been initialized!" );
}

} // end namespace AbstractLinAlgPack
