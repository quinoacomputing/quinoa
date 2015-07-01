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

#include "AbstractLinAlgPack_MatrixPermAggr.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_Permutation.hpp"
#include "AbstractLinAlgPack_PermutationOut.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / initializers

MatrixPermAggr::MatrixPermAggr()
{} // Nothing to explicitly initialize

MatrixPermAggr::MatrixPermAggr(
  const mat_ptr_t      &mat_orig
  ,const perm_ptr_t    &row_perm
  ,const perm_ptr_t    &col_perm
  ,const mat_ptr_t     &mat_perm
  )
{
  this->initialize(mat_orig,row_perm,col_perm,mat_perm);
}

void MatrixPermAggr::initialize(
  const mat_ptr_t      &mat_orig
  ,const perm_ptr_t    &row_perm
  ,const perm_ptr_t    &col_perm
  ,const mat_ptr_t     &mat_perm
  )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    mat_orig.get() == NULL, std::invalid_argument
    ,"MatrixPermAggr::initialize(...): Error!" );
#endif
#ifdef ABSTRACTLINALGPACK_ASSERT_COMPATIBILITY
  bool is_compatible = false;
  if(row_perm.get()) {
    is_compatible = mat_orig->space_cols().is_compatible(row_perm->space());
    TEUCHOS_TEST_FOR_EXCEPTION(
      !is_compatible, VectorSpace::IncompatibleVectorSpaces
      ,"MatrixPermAggr::initialize(...): Error, " 
      "mat_orig->space_cols().is_compatible(row_perm->space()) == false" );
  }
  if(col_perm.get()) {
    is_compatible = mat_orig->space_rows().is_compatible(col_perm->space());
    TEUCHOS_TEST_FOR_EXCEPTION(
      !is_compatible, VectorSpace::IncompatibleVectorSpaces
      ,"MatrixPermAggr::initialize(...): Error, " 
      "mat_orig->space_rows().is_compatible(col_perm->space()) == false" );
  }
#endif
  mat_orig_   = mat_orig;
  row_perm_   = row_perm;
  col_perm_   = col_perm;
  mat_perm_   = mat_perm;
}

void MatrixPermAggr::set_uninitialized()
{
  namespace rcp = MemMngPack;
  mat_orig_   = Teuchos::null;
  row_perm_   = Teuchos::null;
  col_perm_   = Teuchos::null;
  mat_perm_   = Teuchos::null;
}

// Overridden from MatrixBase

size_type MatrixPermAggr::rows() const
{
  return mat_orig_.get() ? mat_orig_->rows() : 0;
}

size_type MatrixPermAggr::cols() const
{
  return mat_orig_.get() ? mat_orig_->cols() : 0;
}

size_type MatrixPermAggr::nz() const
{
  return mat_orig_.get() ? mat_orig_->nz() : 0;
}

// Overridden from MatrixOp

const VectorSpace& MatrixPermAggr::space_cols() const
{
  return mat_orig_->space_cols();
}

const VectorSpace& MatrixPermAggr::space_rows() const
{
  return mat_orig_->space_rows();
}

MatrixOp::mat_ptr_t
MatrixPermAggr::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
  if(mat_perm_.get())
    return mat_perm_->sub_view(row_rng,col_rng);
  if(!row_perm_.get() && !col_perm_.get())
    return mat_orig_->sub_view(row_rng,col_rng);
  return MatrixOp::sub_view(row_rng,col_rng); // ToDo: Speicalized!
}

MatrixOp& MatrixPermAggr::operator=(const MatrixOp& M)
{
  using Teuchos::dyn_cast;
  const MatrixPermAggr
    Mp = dyn_cast<const MatrixPermAggr>(M);
  if( this == &Mp )
    return *this; // Assignment to self
  // Shallow copy is okay as long as client is careful!
  mat_orig_  = Mp.mat_orig_;
  row_perm_  = Mp.row_perm_;
  col_perm_  = Mp.col_perm_;
  mat_perm_  = Mp.mat_perm_;
  return *this;
}

std::ostream& MatrixPermAggr::output(std::ostream& out) const
{
  out << "Matrix with permuted view:\n";
  out << "mat_orig =\n" << *mat_orig_;
  out << "row_perm =";
  if( row_perm_.get() )
    out << "\n" << *row_perm_;
  else
    out << " NULL\n";
  out << "col_perm =";
  if( col_perm_.get() )
    out << "\n" << *col_perm_;
  else
    out << " NULL\n";
  out << "mat_perm =";
  if( mat_perm_.get() )
    out << "\n" << *mat_perm_;
  else
    out << " NULL\n";
  return out;
}

bool MatrixPermAggr::Mp_StM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Mp_StM(mwo_lhs,alpha,*mat_orig_,trans_rhs);
    return true;
  }
  AbstractLinAlgPack::Mp_StM(mwo_lhs,alpha,*mat_orig_,trans_rhs); // ToDo: Specialize!
  return true;
}

bool MatrixPermAggr::Mp_StMtP(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp M_trans
  ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Mp_StMtP(mwo_lhs,alpha,*mat_orig_,M_trans,P_rhs,P_rhs_trans);
    return true;
  }
  AbstractLinAlgPack::Mp_StMtP(mwo_lhs,alpha,*mat_orig_,M_trans,P_rhs,P_rhs_trans); // ToDo: Specialize!
  return true;
}

bool MatrixPermAggr::Mp_StPtM(
  MatrixOp* mwo_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  , BLAS_Cpp::Transp M_trans
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Mp_StPtM(mwo_lhs,alpha,P_rhs,P_rhs_trans,*mat_orig_,M_trans);
    return true;
  }
  AbstractLinAlgPack::Mp_StPtM(mwo_lhs,alpha,P_rhs,P_rhs_trans,*mat_orig_,M_trans); // ToDo: Specialize!
  return true;
}

bool MatrixPermAggr::Mp_StPtMtP(
  MatrixOp* mwo_lhs, value_type alpha
  ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  ,BLAS_Cpp::Transp M_trans
  ,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Mp_StPtMtP(mwo_lhs,alpha,P_rhs1,P_rhs1_trans,*mat_orig_,M_trans,P_rhs2,P_rhs2_trans);
    return true;
  }
  AbstractLinAlgPack::Mp_StPtMtP(mwo_lhs,alpha,P_rhs1,P_rhs1_trans,*mat_orig_,M_trans,P_rhs2,P_rhs2_trans); // ToDo: Specialize!
  return true;
}

void MatrixPermAggr::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  ,const Vector& x, value_type b
  ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using LinAlgOpPack::V_MtV;

  if(mat_perm_.get()) {
    AbstractLinAlgPack::Vp_StMtV(y,a,*mat_perm_,M_trans,x,b);
    return;
  }
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Vp_StMtV(y,a,*mat_orig_,M_trans,x,b);
    return;
  }

  //
  // y = a*op(Pr*M*Pc')*x + b*y
  //
  // =>
  //
  // y = a*P1*op(M)*P2'*x + b*y
  //
  // =>
  //
  // ta = P2'*x
  // tb = op(M)*ta
  // tc = P1*tb
  // y = a*tc + b*y
  //
  const Permutation
    *P1 = ( M_trans == no_trans ? row_perm_.get() : col_perm_.get() ),
    *P2 = ( M_trans == no_trans ? col_perm_.get() : row_perm_.get() );
  VectorSpace::vec_mut_ptr_t  ta, tb, tc;
  // ta = P2'*x
  if( P2 && !P2->is_identity() )
    P2->permute( trans, x, (ta = P2->space().create_member()).get() );
  else
    *(tb = ( M_trans == no_trans ? mat_orig_->space_rows() : mat_orig_->space_cols() ).create_member() ) = x;
  // tb = op(M)*ta
  V_MtV(
    (tb = ( M_trans == no_trans ? mat_orig_->space_cols() : mat_orig_->space_rows() ).create_member() ).get()
    ,*mat_orig_, M_trans, *ta
    );
  // tc = P1*tb
  if( P1 && !P1->is_identity() )
    P1->permute( no_trans, *tb, (tc = P1->space().create_member()).get() );
  else
    tc = tb->clone();
  // y = b*y + a*tc
  AbstractLinAlgPack::Vt_S( y, b );
  AbstractLinAlgPack::Vp_StV( y, a, *tc );
}

void MatrixPermAggr::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp M_trans
  , const SpVectorSlice& x, value_type b) const
{
  if(mat_perm_.get()) {
    AbstractLinAlgPack::Vp_StMtV(y,a,*mat_perm_,M_trans,x,b);
    return;
  }
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Vp_StMtV(y,a,*mat_orig_,M_trans,x,b);
    return;
  }
  MatrixOp::Vp_StMtV(y,a,M_trans,x,b);
}

void MatrixPermAggr::Vp_StPtMtV(
  VectorMutable* vs_lhs, value_type alpha
  ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  ,BLAS_Cpp::Transp M_rhs2_trans
  ,const Vector& v_rhs3, value_type beta
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,*mat_orig_,M_rhs2_trans,v_rhs3,beta);
    return;
  }
  MatrixOp::Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,v_rhs3,beta);
}

void MatrixPermAggr::Vp_StPtMtV(
  VectorMutable* vs_lhs, value_type alpha
  ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  ,BLAS_Cpp::Transp M_rhs2_trans
  ,const SpVectorSlice& sv_rhs3, value_type beta
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,*mat_orig_,M_rhs2_trans,sv_rhs3,beta);
    return;
  }
  MatrixOp::Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,sv_rhs3,beta);
}

value_type MatrixPermAggr::transVtMtV(
  const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
  ,const Vector& v_rhs3
  ) const
{
  if(!row_perm_.get() && !col_perm_.get())
    return AbstractLinAlgPack::transVtMtV(v_rhs1,*mat_orig_,trans_rhs2,v_rhs3);
  return MatrixOp::transVtMtV(v_rhs1,trans_rhs2,v_rhs3);
}

value_type MatrixPermAggr::transVtMtV(
  const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
  ,const SpVectorSlice& sv_rhs3
  ) const
{
  if(!row_perm_.get() && !col_perm_.get())
    return AbstractLinAlgPack::transVtMtV(sv_rhs1,*mat_orig_,trans_rhs2,sv_rhs3);
  return MatrixOp::transVtMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

void MatrixPermAggr::syr2k(
  BLAS_Cpp::Transp M_trans, value_type alpha
  ,const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
  ,const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
  ,value_type beta, MatrixSymOp* symwo_lhs
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::syr2k(*mat_orig_,M_trans,alpha,P1,P1_trans,P2,P2_trans,beta,symwo_lhs);
    return;
  }
  MatrixOp::syr2k(M_trans,alpha,P1,P1_trans,P2,P2_trans,beta,symwo_lhs);
}

bool MatrixPermAggr::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1, const MatrixOp& mwo_rhs2
  ,BLAS_Cpp::Transp trans_rhs2, value_type beta
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Mp_StMtM(mwo_lhs,alpha,*mat_orig_,trans_rhs1,mwo_rhs2,trans_rhs2,beta);
    return true;
  }
  return MatrixOp::Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta);
}

bool MatrixPermAggr::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2, value_type beta
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::Mp_StMtM(mwo_lhs,alpha,mwo_rhs1,trans_rhs1,*mat_orig_,trans_rhs2,beta);
    return true;
  }
  return MatrixOp::Mp_StMtM(mwo_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2,beta);
}

bool MatrixPermAggr::syrk(
  BLAS_Cpp::Transp M_trans, value_type alpha
  ,value_type beta, MatrixSymOp* sym_lhs
  ) const
{
  if(!row_perm_.get() && !col_perm_.get()) {
    AbstractLinAlgPack::syrk(*mat_orig_,M_trans,alpha,beta,sym_lhs);
    return true;
  }
  return MatrixOp::syrk(M_trans,alpha,beta,sym_lhs);
}

} // end namespace AbstractLinAlgPack
