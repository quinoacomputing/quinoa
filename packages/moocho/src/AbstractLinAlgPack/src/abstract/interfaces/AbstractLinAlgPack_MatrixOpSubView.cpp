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

#include <typeinfo>
#include <stdexcept>

#include "AbstractLinAlgPack_MatrixOpSubView.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorView.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

namespace AbstractLinAlgPack {

MatrixOpSubView::MatrixOpSubView(
  const mat_ptr_t& M_full
  ,const Range1D& rng_rows
  ,const Range1D& rng_cols
  ,BLAS_Cpp::Transp M_trans
  )
{
  this->initialize(M_full,rng_rows,rng_cols,M_trans);
}
    
void MatrixOpSubView::initialize(
  const mat_ptr_t& M_full
  ,const Range1D& rng_rows_in
  ,const Range1D& rng_cols_in
  ,BLAS_Cpp::Transp M_trans
  )
{
  namespace rcp = MemMngPack;

  if( M_full.get() ) {
    const index_type
      M_rows = M_full->rows(),
      M_cols = M_full->cols();
    const Range1D
      rng_rows = RangePack::full_range(rng_rows_in,1,M_rows),
      rng_cols = RangePack::full_range(rng_cols_in,1,M_cols);
    TEUCHOS_TEST_FOR_EXCEPTION(
      rng_rows.ubound() > M_rows, std::invalid_argument
      ,"MatrixOpSubView::initialize(...): Error, "
      "rng_rows = ["<<rng_rows.lbound()<<","<<rng_rows.ubound()<<"] is of range of "
      "[1,M_full->rows()] = [1,"<<M_rows<<"]" );
    TEUCHOS_TEST_FOR_EXCEPTION(
      rng_cols.ubound() > M_cols, std::invalid_argument
      ,"MatrixOpSubView::initialize(...): Error, "
      "rng_cols = ["<<rng_cols.lbound()<<","<<rng_cols.ubound()<<"] is of range of "
      "[1,M_full->cols()] = [1,"<<M_cols<<"]" );
    M_full_     = M_full;
    rng_rows_   = rng_rows;
    rng_cols_   = rng_cols;
    M_trans_    = M_trans;
    space_cols_ = ( M_trans == BLAS_Cpp::no_trans
            ? M_full->space_cols().sub_space(rng_rows)->clone()
            : M_full->space_rows().sub_space(rng_cols)->clone() );
    space_rows_ = ( M_trans == BLAS_Cpp::no_trans
            ? M_full->space_rows().sub_space(rng_cols)->clone()
            : M_full->space_cols().sub_space(rng_rows)->clone() );
  }
  else {
    M_full_     = Teuchos::null;
    rng_rows_   = Range1D::Invalid;
    rng_cols_   = Range1D::Invalid;
    M_trans_    = BLAS_Cpp::no_trans;
    space_cols_ = Teuchos::null;
    space_rows_ = Teuchos::null;
  }
}

// overridden from MatrixBase

size_type MatrixOpSubView::rows() const
{
  return ( M_full_.get() 
       ? BLAS_Cpp::rows( rng_rows_.size(), rng_cols_.size(), M_trans_ )
       : 0 );
}

size_type MatrixOpSubView::cols() const
{
  return ( M_full_.get() 
       ? BLAS_Cpp::cols( rng_rows_.size(), rng_cols_.size(), M_trans_ )
       : 0 );
}

size_type MatrixOpSubView::nz() const
{
  return ( M_full_.get()
       ? ( rng_rows_.full_range() && rng_cols_.full_range()
         ? M_full_->nz()
         : MatrixBase::nz() )
       : 0 );
}

// Overridden form MatrixOp

const VectorSpace& MatrixOpSubView::space_cols() const
{
  assert_initialized();
  return *space_cols_;
}

const VectorSpace& MatrixOpSubView::space_rows() const
{
  assert_initialized();
  return *space_rows_;
}

MatrixOp::mat_ptr_t
MatrixOpSubView::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
  assert_initialized();
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
  return Teuchos::null;
}

void MatrixOpSubView::zero_out()
{
  assert_initialized();
  if( rng_rows_.full_range() && rng_cols_.full_range() ) {
    M_full_->zero_out();
    return;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error, "MatrixOpSubView::zero_out(): "
    "Error, this method can not be implemented with a sub-view" );
}

void MatrixOpSubView::Mt_S( value_type alpha )
{
  assert_initialized();
  if( rng_rows_.full_range() && rng_cols_.full_range() ) {
    M_full_->Mt_S(alpha);
    return;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error, "MatrixOpSubView::Mt_S(alpha): "
    "Error, this method can not be implemented with a sub-view" );
}

MatrixOp& MatrixOpSubView::operator=(const MatrixOp& M)
{
  assert_initialized();
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
  return *this;
}

std::ostream& MatrixOpSubView::output(std::ostream& out) const
{
  assert_initialized();
  return MatrixOp::output(out); // ToDo: Specialize if needed?
}

// Level-1 BLAS

// rhs matrix argument

bool MatrixOpSubView::Mp_StM(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs) const
{
  assert_initialized();
  return MatrixOp::Mp_StM(m_lhs,alpha,trans_rhs); // ToDo: Specialize?
}

bool MatrixOpSubView::Mp_StMtP(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  ) const
{
  assert_initialized();
  return MatrixOp::Mp_StMtP(m_lhs,alpha,M_trans,P_rhs,P_rhs_trans); // ToDo: Specialize?
}

bool MatrixOpSubView::Mp_StPtM(
  MatrixOp* m_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  , BLAS_Cpp::Transp M_trans
  ) const
{
  assert_initialized();
  return MatrixOp::Mp_StPtM(m_lhs,alpha,P_rhs,P_rhs_trans,M_trans); // ToDo: Specialize?
}

bool MatrixOpSubView::Mp_StPtMtP(
  MatrixOp* m_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  ) const
{
  assert_initialized();
  return MatrixOp::Mp_StPtMtP(
    m_lhs,alpha,P_rhs1,P_rhs1_trans,M_trans,P_rhs2,P_rhs2_trans); // ToDo: Specialize?
}

// lhs matrix argument

bool MatrixOpSubView::Mp_StM(
  value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
  assert_initialized();
  return MatrixOp::Mp_StM(alpha,M_rhs,trans_rhs); // ToDo: Specialize?
}

bool MatrixOpSubView::Mp_StMtP(
  value_type alpha
  ,const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
  ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  )
{
  assert_initialized();
  return MatrixOp::Mp_StMtP(alpha,M_rhs,M_trans,P_rhs,P_rhs_trans); // ToDo: Specialize?
}

bool MatrixOpSubView::Mp_StPtM(
  value_type alpha
  ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  ,const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
  )
{
  assert_initialized();
  return MatrixOp::Mp_StPtM(
    alpha,P_rhs,P_rhs_trans,M_rhs,M_trans); // ToDo: Specialize?
}

bool MatrixOpSubView::Mp_StPtMtP(
  value_type alpha
  ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  ,const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
  ,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  )
{
  assert_initialized();
  return MatrixOp::Mp_StPtMtP(
    alpha,P_rhs1,P_rhs1_trans,M_rhs,M_trans,P_rhs2,P_rhs2_trans); // ToDo: Specialize?
}

// Level-2 BLAS

void MatrixOpSubView::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans_in
  , const Vector& x, value_type b
  ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;

  assert_initialized();

  BLAS_Cpp::Transp
    M_trans_trans = ( M_trans_==no_trans ? M_trans_in : BLAS_Cpp::trans_not(M_trans_in) );

  // ToDo: Assert compatibility!

  if( rng_rows_.full_range() && rng_cols_.full_range() ) {
    AbstractLinAlgPack::Vp_StMtV(  // The matrix is just transposed
      y, a
      ,*M_full_, M_trans_trans
      ,x, b );
    return;
  }
  // y = b*y
  Vt_S( y, b );
  //
  // xt1                      = 0.0
  // xt3 = xt(op_op_rng_cols) = x
  // xt3                      = 0.0
  //
  // [ yt1 ]                        [ xt1 ]
  // [ yt2 ] = a * op(op(M_full)) * [ xt3 ]
  // [ yt3 ]                        [ xt3 ]
  //
  // =>
  //
  // y += yt2 = yt(op_op_rng_rows)
  //
  const Range1D
    op_op_rng_rows = ( M_trans_trans == no_trans ? rng_rows_ : rng_cols_ ),
    op_op_rng_cols = ( M_trans_trans == no_trans ? rng_cols_ : rng_rows_ );
  VectorSpace::vec_mut_ptr_t
    xt = ( M_trans_trans == no_trans ? M_full_->space_rows() : M_full_->space_cols() ).create_member(),
    yt = ( M_trans_trans == no_trans ? M_full_->space_cols() : M_full_->space_rows() ).create_member();
  *xt = 0.0;
  *xt->sub_view(op_op_rng_cols) = x;
    LinAlgOpPack::V_StMtV( yt.get(), a, *M_full_, M_trans_trans, *xt );
  LinAlgOpPack::Vp_V( y, *yt->sub_view(op_op_rng_rows) );
}

void MatrixOpSubView::Vp_StMtV(
  VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2, value_type beta) const
{
  assert_initialized();
  MatrixOp::Vp_StMtV(v_lhs,alpha,trans_rhs1,sv_rhs2,beta); // ToDo: Specialize?
}

void MatrixOpSubView::Vp_StPtMtV(
  VectorMutable* v_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const Vector& v_rhs3, value_type beta) const
{
  assert_initialized();
  MatrixOp::Vp_StPtMtV(
    v_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,v_rhs3,beta); // ToDo: Specialize?
}

void MatrixOpSubView::Vp_StPtMtV(
  VectorMutable* v_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const SpVectorSlice& sv_rhs3, value_type beta) const
{
  assert_initialized();
  MatrixOp::Vp_StPtMtV(
    v_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,sv_rhs3,beta); // ToDo: Specialize?
}

value_type MatrixOpSubView::transVtMtV(
  const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const Vector& v_rhs3) const
{
  assert_initialized();
  return MatrixOp::transVtMtV(v_rhs1,trans_rhs2,v_rhs3); // ToDo: Specialize?
}

value_type MatrixOpSubView::transVtMtV(
  const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const SpVectorSlice& sv_rhs3) const
{
  assert_initialized();
  return MatrixOp::transVtMtV(sv_rhs1,trans_rhs2,sv_rhs3); // ToDo: Specialize?
}

void MatrixOpSubView::syr2k(
  BLAS_Cpp::Transp M_trans, value_type alpha
  , const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
  , const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
  , value_type beta, MatrixSymOp* sym_lhs ) const
{
  assert_initialized();
  MatrixOp::syr2k(
    M_trans,alpha,P1,P1_trans,P2,P2_trans,beta,sym_lhs); // ToDo: Specialize?
}

// Level-3 BLAS

bool MatrixOpSubView::Mp_StMtM(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOp& mwo_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
  assert_initialized();
  return MatrixOp::Mp_StMtM(
    m_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta); // ToDo: Specialize?
}

bool MatrixOpSubView::Mp_StMtM(
  MatrixOp* m_lhs, value_type alpha
  , const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  , BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
  return MatrixOp::Mp_StMtM(
    m_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2,beta); // ToDo: Specialize?
}

bool MatrixOpSubView::Mp_StMtM(
  value_type alpha
  ,const MatrixOp& mvw_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
  ,value_type beta )
{
  assert_initialized();
  return MatrixOp::Mp_StMtM(
    alpha,mvw_rhs1,trans_rhs1,mwo_rhs2,trans_rhs2,beta); // ToDo: Specialize?
}

bool MatrixOpSubView::syrk(
  BLAS_Cpp::Transp M_trans, value_type alpha
  , value_type beta, MatrixSymOp* sym_lhs ) const
{
  assert_initialized();
  return MatrixOp::syrk(M_trans,alpha,beta,sym_lhs); // ToDo: Specialize?
}

// private

void MatrixOpSubView::assert_initialized() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    M_full_.get() == NULL, std::logic_error
    ,"Error, the MatrixOpSubView object has not been initialize!" );
}

} // end namespace AbstractLinAlgPack
