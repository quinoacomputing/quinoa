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
#include <math.h>

#include <typeinfo>
#include <stdexcept>

#include "AbstractLinAlgPack_MatrixSymOp.hpp"
#include "AbstractLinAlgPack_MatrixOpSubView.hpp"
#include "AbstractLinAlgPack_MatrixPermAggr.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_Permutation.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorView.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"

namespace AbstractLinAlgPack {

void MatrixOp::zero_out()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error, "MatrixOp::zero_out(): "
    "Error, this method as not been defined by the subclass \'"
    <<typeName(*this)<<"\'" );
}

void MatrixOp::Mt_S(value_type alpha)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error, "MatrixOp::Mt_S(): "
    "Error, this method as not been defined by the subclass \'"
    <<typeName(*this)<<"\'" );
}

MatrixOp& MatrixOp::operator=(const MatrixOp& M)
{
  const bool assign_to_self = dynamic_cast<const void*>(this) == dynamic_cast<const void*>(&M);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !assign_to_self, std::logic_error
    ,"MatrixOp::operator=(M) : Error, this is not assignment "
    "to self and this method is not overridden for the subclass \'"
    << typeName(*this) << "\'" );
  return *this; // assignment to self
}

std::ostream& MatrixOp::output(std::ostream& out_arg) const
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream(Teuchos::rcp(&out_arg,false));
  Teuchos::OSTab tab(out);
  const size_type m = this->rows(), n = this->cols();
  VectorSpace::vec_mut_ptr_t
    row_vec = space_rows().create_member(); // dim() == n
  *out << m << " " << n << std::endl;
  for( size_type i = 1; i <= m; ++i ) {
    LinAlgOpPack::V_StMtV( &*row_vec, 1.0, *this, BLAS_Cpp::trans, EtaVector(i,m)() );
    row_vec->output(*out,false,true);
  }
  return out_arg;
}

// Clone

MatrixOp::mat_mut_ptr_t
MatrixOp::clone()
{
  return Teuchos::null;
}

MatrixOp::mat_ptr_t
MatrixOp::clone() const
{
  return Teuchos::null;
}

// Norms

const MatrixOp::MatNorm
MatrixOp::calc_norm(
  EMatNormType  requested_norm_type
  ,bool         allow_replacement
  ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using LinAlgOpPack::V_MtV;
  const VectorSpace
    &space_cols = this->space_cols(),
    &space_rows = this->space_rows();
  const index_type
    num_rows = space_cols.dim(),
    num_cols = space_rows.dim();
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(requested_norm_type == MAT_NORM_1 || requested_norm_type == MAT_NORM_INF), MethodNotImplemented
    ,"MatrixOp::calc_norm(...): Error, This default implemenation can only "
    "compute the one norm or the infinity norm!"
    );
  //
  // Here we implement Algorithm 2.5 in "Applied Numerical Linear Algebra", Demmel (1997)
  // using the momenclature in the text.
  //
  const MatrixOp
    &B = *this;
  bool
    do_trans = requested_norm_type == MAT_NORM_INF;
  VectorSpace::vec_mut_ptr_t
    x    = (!do_trans ? space_rows : space_cols).create_member(1.0/(!do_trans ? num_cols : num_rows)),
    w    = (!do_trans ? space_cols : space_rows).create_member(),
    zeta = (!do_trans ? space_cols : space_rows).create_member(),
    z    = (!do_trans ? space_rows : space_cols).create_member();
  const index_type max_iter = 5; // Recommended by Highman 1988, (see Demmel's reference)
  value_type w_nrm = 0.0;
  for( index_type k = 0; k <= max_iter; ++k ) {
    V_MtV( w.get(), B, !do_trans ? no_trans : trans, *x );    // w = B*x
    sign( *w, zeta.get() );                                   // zeta = sign(w)
    V_MtV( z.get(), B, !do_trans ? trans : no_trans, *zeta ); // z = B'*zeta
    value_type  z_j = 0.0;                                    // max |z(j)| = ||z||inf
    index_type  j   = 0;
    max_abs_ele( *z, &z_j, &j );
    const value_type zTx = dot(*z,*x);                        // z'*x
    w_nrm = w->norm_1();                                      // ||w||1
    if( ::fabs(z_j) <= zTx ) {                                // Update
      break;
    }
    else {
      *x = 0.0;
      x->set_ele(j,1.0);
    }
  }
  return MatNorm(w_nrm,requested_norm_type);
}

// Subview

MatrixOp::mat_ptr_t
MatrixOp::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
  namespace rcp = MemMngPack;

  if( 
    ( ( row_rng.lbound() == 1 && row_rng.ubound() == this->rows() )
      || row_rng.full_range() )
    &&
    ( ( col_rng.lbound() == 1 && col_rng.ubound() == this->cols() )
      || row_rng.full_range() )
    ) 
  {
    return Teuchos::rcp(this,false); // don't clean up memory
  }
  return Teuchos::rcp(
    new MatrixOpSubView(
      Teuchos::rcp(const_cast<MatrixOp*>(this),false) // don't clean up memory
      ,row_rng,col_rng ) );
}

// Permuted views

MatrixOp::mat_ptr_t
MatrixOp::perm_view(
  const Permutation          *P_row
  ,const index_type          row_part[]
  ,int                       num_row_part
  ,const Permutation         *P_col
  ,const index_type          col_part[]
  ,int                       num_col_part
  ) const
{
  namespace rcp = MemMngPack;
  return Teuchos::rcp(
    new MatrixPermAggr(
      Teuchos::rcp(this,false)
      ,Teuchos::rcp(P_row,false)
      ,Teuchos::rcp(P_col,false)
      ,Teuchos::null
      ) );
}

MatrixOp::mat_ptr_t
MatrixOp::perm_view_update(
  const Permutation          *P_row
  ,const index_type          row_part[]
  ,int                       num_row_part
  ,const Permutation         *P_col
  ,const index_type          col_part[]
  ,int                       num_col_part
  ,const mat_ptr_t           &perm_view
  ) const
{
  return this->perm_view(
    P_row,row_part,num_row_part
    ,P_col,col_part,num_col_part );
}

// Level-1 BLAS

bool MatrixOp::Mp_StM(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs) const
{
  return false;
}

bool MatrixOp::Mp_StM(
  value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
  return false;
}

bool MatrixOp::Mp_StMtP(
  MatrixOp* m_lhs, value_type alpha
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  ) const
{
  return false;
}

bool MatrixOp::Mp_StMtP(
  value_type alpha
  ,const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
  ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  )
{
  return false;
}

bool MatrixOp::Mp_StPtM(
  MatrixOp* m_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  , BLAS_Cpp::Transp M_trans
  ) const
{
  return false;
}

bool MatrixOp::Mp_StPtM(
  value_type alpha
  ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  ,const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
  )
{
  return false;
}

bool MatrixOp::Mp_StPtMtP(
  MatrixOp* m_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  ) const
{
  return false;
}

bool MatrixOp::Mp_StPtMtP(
  value_type alpha
  ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  ,const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
  ,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  )
{
  return false;
}

// Level-2 BLAS

void MatrixOp::Vp_StMtV(
  VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2, value_type beta) const
{
  Vp_MtV_assert_compatibility(v_lhs,*this,trans_rhs1,sv_rhs2 );
  if( sv_rhs2.nz() ) {
    VectorSpace::vec_mut_ptr_t
      v_rhs2 = (trans_rhs1 == BLAS_Cpp::no_trans
            ? this->space_rows()
            : this->space_cols()
        ).create_member();
    v_rhs2->set_sub_vector(sub_vec_view(sv_rhs2));
    this->Vp_StMtV(v_lhs,alpha,trans_rhs1,*v_rhs2,beta);
  }
  else {
    Vt_S( v_lhs, beta );
  }
}

void MatrixOp::Vp_StPtMtV(
  VectorMutable* y, value_type a
  ,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  ,BLAS_Cpp::Transp M_trans
  ,const Vector& x, value_type b
  ) const
{
  VectorSpace::vec_mut_ptr_t
    t = ( M_trans == BLAS_Cpp::no_trans
        ? this->space_cols()
        : this->space_rows() ).create_member();
  LinAlgOpPack::V_MtV( t.get(), *this, M_trans, x );
  AbstractLinAlgPack::Vp_StMtV( y, a, P, P_trans, *t, b );
}

void MatrixOp::Vp_StPtMtV(
  VectorMutable* y, value_type a
  ,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  ,BLAS_Cpp::Transp M_trans
  ,const SpVectorSlice& x, value_type b
  ) const
{
  VectorSpace::vec_mut_ptr_t
    t = ( M_trans == BLAS_Cpp::no_trans
        ? this->space_cols()
        : this->space_rows() ).create_member();
  LinAlgOpPack::V_MtV( t.get(), *this, M_trans, x );
  AbstractLinAlgPack::Vp_StMtV( y, a, P, P_trans, *t, b );
}

value_type MatrixOp::transVtMtV(
  const Vector& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const Vector& vs_rhs3) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
  return 0.0;
}

value_type MatrixOp::transVtMtV(
  const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const SpVectorSlice& sv_rhs3) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
  return 0.0;
}

void MatrixOp::syr2k(
  BLAS_Cpp::Transp M_trans, value_type alpha
  , const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
  , const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
  , value_type beta, MatrixSymOp* sym_lhs ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
}

// Level-3 BLAS

bool MatrixOp::Mp_StMtM(
  MatrixOp* C, value_type a
  ,BLAS_Cpp::Transp A_trans, const MatrixOp& B
  ,BLAS_Cpp::Transp B_trans, value_type b) const
{
  return false;
}

bool MatrixOp::Mp_StMtM(
  MatrixOp* m_lhs, value_type alpha
  ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
  return false;
}

bool MatrixOp::Mp_StMtM(
  value_type alpha
  ,const MatrixOp& mvw_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
  ,value_type beta )
{
  return false;
}

bool MatrixOp::syrk(
  BLAS_Cpp::Transp   M_trans
  ,value_type        alpha
  ,value_type        beta
  ,MatrixSymOp   *sym_lhs
  ) const
{
  return false;
}

bool MatrixOp::syrk(
  const MatrixOp  &mwo_rhs
  ,BLAS_Cpp::Transp   M_trans
  ,value_type         alpha
  ,value_type         beta
  )
{
  return false;
}

} // end namespace AbstractLinAlgPack

// Non-member functions

// level-1 BLAS

void AbstractLinAlgPack::Mt_S( MatrixOp* m_lhs, value_type alpha )
{
  if(alpha == 0.0)
    m_lhs->zero_out();
  else if( alpha != 1.0 )
    m_lhs->Mt_S(alpha);
}

void AbstractLinAlgPack::Mp_StM(
  MatrixOp* mwo_lhs, value_type alpha, const MatrixOp& M_rhs
  , BLAS_Cpp::Transp trans_rhs)
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;

  // Give the rhs argument a chance to implement the operation
  if(M_rhs.Mp_StM(mwo_lhs,alpha,trans_rhs))
    return;

  // Give the lhs argument a change to implement the operation
  if(mwo_lhs->Mp_StM(alpha,M_rhs,trans_rhs))
    return;

  // We must try to implement the method
  MultiVectorMutable
    *m_mut_lhs = dynamic_cast<MultiVectorMutable*>(mwo_lhs);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !m_mut_lhs || !(m_mut_lhs->access_by() & MultiVector::COL_ACCESS)
    ,MatrixOp::MethodNotImplemented
    ,"MatrixOp::Mp_StM(...) : Error, mwo_lhs of type \'"
    << typeName(*mwo_lhs) << "\' could not implement the operation "
    "and does not support the "
    "\'MultiVectorMutable\' interface.  Furthermore, "
    "the rhs matix argument of type \'" << typeName(*mwo_lhs)
    << "\' could not implement the operation!" );
    
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    !mwo_lhs->space_rows().is_compatible(
      trans_rhs == no_trans ? M_rhs.space_rows() : M_rhs.space_cols() )
    || !mwo_lhs->space_cols().is_compatible(
      trans_rhs == no_trans ? M_rhs.space_cols() : M_rhs.space_rows() )
    , MatrixOp::IncompatibleMatrices
    ,"MatrixOp::Mp_StM(mwo_lhs,...): Error, mwo_lhs of type \'"<<typeName(*mwo_lhs)<<"\' "
    <<"is not compatible with M_rhs of type \'"<<typeName(M_rhs)<<"\'" );
#endif

  const size_type
    rows = BLAS_Cpp::rows( mwo_lhs->rows(), mwo_lhs->cols(), trans_rhs ),
    cols = BLAS_Cpp::cols( mwo_lhs->rows(), mwo_lhs->cols(), trans_rhs );
  for( size_type j = 1; j <= cols; ++j )
    AbstractLinAlgPack::Vp_StMtV( m_mut_lhs->col(j).get(), alpha, M_rhs, trans_rhs, EtaVector(j,cols)() );
  // ToDo: consider row access?
}

void AbstractLinAlgPack::Mp_StMtP(
  MatrixOp* mwo_lhs, value_type alpha
  , const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  )
{

  // Give the rhs argument a chance to implement the operation
  if(M_rhs.Mp_StMtP(mwo_lhs,alpha,M_trans,P_rhs,P_rhs_trans))
    return;

  // Give the lhs argument a change to implement the operation
  if(mwo_lhs->Mp_StMtP(alpha,M_rhs,M_trans,P_rhs,P_rhs_trans))
    return;

  // We must try to implement the method
  MultiVectorMutable
    *m_mut_lhs = dynamic_cast<MultiVectorMutable*>(mwo_lhs);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !m_mut_lhs, MatrixOp::MethodNotImplemented
    ,"MatrixOp::Mp_StMtP(...) : Error, mwo_lhs of type \'"
    << typeName(*mwo_lhs) << "\' does not support the "
    "\'MultiVectorMutable\' interface!" );

  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
}

void AbstractLinAlgPack::Mp_StPtM(
  MatrixOp* mwo_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  , const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
  )
{

  // Give the rhs argument a chance to implement the operation
  if(M_rhs.Mp_StPtM(mwo_lhs,alpha,P_rhs,P_rhs_trans,M_trans))
    return;

  // Give the lhs argument a change to implement the operation
  if(mwo_lhs->Mp_StPtM(alpha,P_rhs,P_rhs_trans,M_rhs,M_trans))
    return;

  // We must try to implement the method
  MultiVectorMutable
    *m_mut_lhs = dynamic_cast<MultiVectorMutable*>(mwo_lhs);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !m_mut_lhs, MatrixOp::MethodNotImplemented
    ,"MatrixOp::Mp_StPtM(...) : Error, mwo_lhs of type \'"
    << typeName(*mwo_lhs) << "\' does not support the "
    "\'MultiVectorMutable\' interface!" );

  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!

}

void AbstractLinAlgPack::Mp_StPtMtP(
  MatrixOp* mwo_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs
  , const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  )
{

  // Give the rhs argument a chance to implement the operation
  if(M_rhs.Mp_StPtMtP(mwo_lhs,alpha,P_rhs1,P_rhs1_trans,trans_rhs,P_rhs2,P_rhs2_trans))
    return;

  // Give the lhs argument a change to implement the operation
  if(mwo_lhs->Mp_StPtMtP(alpha,P_rhs1,P_rhs1_trans,M_rhs,trans_rhs,P_rhs2,P_rhs2_trans))
    return;

  // We must try to implement the method
  MultiVectorMutable
    *m_mut_lhs = dynamic_cast<MultiVectorMutable*>(mwo_lhs);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !m_mut_lhs, MatrixOp::MethodNotImplemented
    ,"MatrixOp::Mp_StPtMtP(...) : Error, mwo_lhs of type \'"
    << typeName(*mwo_lhs) << "\' does not support the "
    "\'MultiVectorMutable\' interface!" );

  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!

}

// level-3 blas

void AbstractLinAlgPack::Mp_StMtM(
  MatrixOp* C, value_type a
  ,const MatrixOp& A, BLAS_Cpp::Transp A_trans
  ,const MatrixOp& B, BLAS_Cpp::Transp B_trans
  ,value_type b
  )
{
  
  // Give A a chance
  if(A.Mp_StMtM(C,a,A_trans,B,B_trans,b))
    return;
  // Give B a chance
  if(B.Mp_StMtM(C,a,A,A_trans,B_trans,b))
    return;
  // Give C a chance
  if(C->Mp_StMtM(a,A,A_trans,B,B_trans,b))
    return;

  //
  // C = b*C + a*op(A)*op(B)
  //
  // We will perform this by column as:
  //
  //   C(:,j) = b*C(:,j) + a*op(A)*op(B)*e(j), for j = 1...C.cols()
  //
  //   by performing:
  //
  //        t = op(B)*e(j)
  //        C(:,j) = b*C(:,j) + a*op(A)*t
  //
  Mp_MtM_assert_compatibility(C,BLAS_Cpp::no_trans,A,A_trans,B,B_trans);
  MultiVectorMutable *Cmv = dynamic_cast<MultiVectorMutable*>(C);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !Cmv || !(Cmv->access_by() & MultiVector::COL_ACCESS)
    ,MatrixOp::MethodNotImplemented
    ,"AbstractLinAlgPack::Mp_StMtM(...) : Error, mwo_lhs of type \'"
    << typeName(*C) << "\' does not support the "
    "\'MultiVectorMutable\' interface or does not support column access!" );
  // ToDo: We could do this by row also!
  VectorSpace::vec_mut_ptr_t
    t = ( B_trans == BLAS_Cpp::no_trans ? B.space_cols() : B.space_rows() ).create_member();
  const index_type
    C_rows = Cmv->rows(),
    C_cols = Cmv->cols();
  for( index_type j = 1; j <= C_cols; ++j ) {
    // t = op(B)*e(j)
    LinAlgOpPack::V_MtV( t.get(), B, B_trans, EtaVector(j,C_cols)	);	
    // C(:,j) = a*op(A)*t + b*C(:,j)
    Vp_StMtV( Cmv->col(j).get(), a, A, A_trans, *t, b );
  }
}

void AbstractLinAlgPack::syrk(
  const MatrixOp  &A
  ,BLAS_Cpp::Transp   A_trans
  ,value_type         a
  ,value_type         b
  ,MatrixSymOp    *B
  )
{
  // Give A a chance
  if(A.syrk(A_trans,a,b,B))
    return;
  // Give B a chance
  if(B->syrk(A,A_trans,a,b))
    return;
  
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, MatrixOp::MethodNotImplemented
    ,"AbstractLinAlgPack::syrk(...) : Error, neither the right-hand-side matrix "
    "argument mwo_rhs of type \'" << typeName(A) << " nore the left-hand-side matrix "
    "argument sym_lhs of type \'" << typeName(*B) << "\' could implement this operation!"
    );

}
