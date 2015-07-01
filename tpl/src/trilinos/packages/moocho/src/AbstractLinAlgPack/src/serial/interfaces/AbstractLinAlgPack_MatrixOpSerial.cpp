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

// ToDo: 7/27/99: Give these default implementations and test them.

#include <typeinfo>

#include "AbstractLinAlgPack_MatrixOpSerial.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_MatrixOpGetGMSMutable.hpp"
#include "AbstractLinAlgPack_MatrixOpGetGMSTri.hpp"
#include "AbstractLinAlgPack_MatrixSymOpGetGMSSymMutable.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::Vp_StMtV;
  using AbstractLinAlgPack::Mp_StM;
}

namespace AbstractLinAlgPack {

// Level-1 BLAS

void MatrixOpSerial::Mp_StM(DMatrixSlice* gms_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs) const
{
  DenseLinAlgPack::Mp_M_assert_sizes( gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
    , rows(), cols(), trans_rhs );
  const size_type
    m = gms_lhs->rows(),
    n = gms_lhs->cols();
  //
  // Use sparse matrix-vector multiplication to perform this operation.
  // C += a * B = a * B * I = [ a*B*e(1), a*B*e(2), ..., a*B*e(m) ]
  //
  SpVector rhs;
  rhs.uninitialized_resize( n, 1, 1 );
  for( size_type j = 1; j <=n; ++j ) {
    rhs.begin()->initialize( j, 1.0 );	// e(j)
    this->Vp_StMtV( &gms_lhs->col(j), alpha, trans_rhs, rhs(), 1.0 );
  }
}

void MatrixOpSerial::Mp_StMtP(DMatrixSlice* C, value_type a
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  ) const 
{
  // C += a * op(M) * op(P)
  TEUCHOS_TEST_FOR_EXCEPT(true);	// Implement this!
}

void MatrixOpSerial::Mp_StPtM(DMatrixSlice* C, value_type a
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp M_trans
  ) const 
{
  // C += a * op(P) * op(M)
  TEUCHOS_TEST_FOR_EXCEPT(true);	// Implement this!
}

void MatrixOpSerial::Mp_StPtMtP( DMatrixSlice* C, value_type a
  , const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
  ) const
{
  // C += a * op(P1) * op(M) * op(P2)
  TEUCHOS_TEST_FOR_EXCEPT(true);	// Implement this!
}

// Level-2 BLAS

void MatrixOpSerial::Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2, value_type beta) const
{
  Vp_MtV_assert_sizes( vs_lhs->dim(), rows(), cols(), trans_rhs1, sv_rhs2.dim() );
  if( !sv_rhs2.nz() ) {
    // vs_lhs = beta * vs_lhs
    if(beta==0.0)      *vs_lhs = 0.0;
    else if(beta!=1.0) DenseLinAlgPack::Vt_S(vs_lhs,beta);
  }
  else {
    // Convert to dense by default.
    if( sv_rhs2.dim() == sv_rhs2.nz() && sv_rhs2.is_sorted() ) {
      const DVectorSlice vs_rhs2 = AbstractLinAlgPack::dense_view(sv_rhs2);
      this->Vp_StMtV( vs_lhs, alpha, trans_rhs1, vs_rhs2, beta );
    }
    else {
      DVector v_rhs2;
      LinAlgOpPack::assign( &v_rhs2, sv_rhs2 );
      this->Vp_StMtV( vs_lhs, alpha, trans_rhs1, v_rhs2(), beta );
    }
  }
}

void MatrixOpSerial::Vp_StPtMtV(DVectorSlice* y, value_type a
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp M_trans
  , const DVectorSlice& x, value_type b) const
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  Workspace<value_type> t_ws(wss,BLAS_Cpp::cols(P.rows(),P.cols(),P_trans));
  DVectorSlice                t(&t_ws[0],t_ws.size());
    LinAlgOpPack::V_StMtV(&t,a,*this,M_trans,x);
  LinAlgOpPack::Vp_MtV( y, P, P_trans, t, b ); 
}

void MatrixOpSerial::Vp_StPtMtV(DVectorSlice* y, value_type a
  , const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  , BLAS_Cpp::Transp M_trans
  , const SpVectorSlice& x, value_type b) const
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  Workspace<value_type> t_ws(wss,BLAS_Cpp::cols(P.rows(),P.cols(),P_trans));
  DVectorSlice                t(&t_ws[0],t_ws.size());
    LinAlgOpPack::V_StMtV(&t,a,*this,M_trans,x);
  LinAlgOpPack::Vp_MtV( y, P, P_trans, t, b ); 
}

value_type MatrixOpSerial::transVtMtV(const DVectorSlice& x1
  , BLAS_Cpp::Transp M_trans, const DVectorSlice& x2) const
{
  DenseLinAlgPack::Vp_MtV_assert_sizes( x1.dim(), rows(), cols(), M_trans, x2.dim() );
  DVector tmp(x1.dim());
  this->Vp_StMtV( &tmp(), 1.0, M_trans, x2, 0.0 );
  return DenseLinAlgPack::dot( x1, tmp() );
}

value_type MatrixOpSerial::transVtMtV(const SpVectorSlice& x1
  , BLAS_Cpp::Transp M_trans, const SpVectorSlice& x2) const
{
  DenseLinAlgPack::Vp_MtV_assert_sizes( x1.dim(), rows(), cols(), M_trans, x2.dim() );
  if( !x1.nz() || !x2.nz() ) {
    return 0.0;
  }
  else {
    if( x1.overlap(x2) == DenseLinAlgPack::SAME_MEM && x1.dim() == x1.nz() && x1.is_sorted()  ) {
      const DVectorSlice x1_d = AbstractLinAlgPack::dense_view(x1);
      return this->transVtMtV( x1_d, M_trans, x1_d );
    }
    DVector tmp(x1.dim());
    this->Vp_StMtV( &tmp(), 1.0, M_trans, x2, 0.0 );
    return AbstractLinAlgPack::dot( x1, tmp() );
  }
}

void MatrixOpSerial::syr2k(
   BLAS_Cpp::Transp M_trans_in, value_type a
  , const GenPermMatrixSlice& P1_in, BLAS_Cpp::Transp P1_trans
  , const GenPermMatrixSlice& P2_in, BLAS_Cpp::Transp P2_trans
  , value_type b, DMatrixSliceSym* S ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::trans_not;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  //
  // S = b * S
  //
  // S += a*op(P1')*op(M)*op(P2) + a*op(P2')*op(M')*op(P1)
  //
  // We will start by renaming P1 and P2 such that op(P1).rows() >= op(P2).rows().
  // This is because we are going to store some temparary vectors and we don't
  // want them to be too big.
  //
  // We will perform the above operation by working with columns of:
  //
  //    op(P1)(:,j(k)) = e(i(k)) <: R^n
  //
  // Then for each column in op(P1) we will perform:
  //
  //
  // for k = 1...P1.nz()
  //
  //              [    .     ]
  //     S += a * [ e(i(k))' ] * op(M)*op(P2) + a * op(P2') * op(M') * [  ...  e(i(k))  ...  ] 
  //              [    .     ]
  //                row j(k)                                                    col j(k)
  //     =>
  //              [  .   ]
  //     S += a * [ y_k' ] + a * [  ...  y_k  ...  ] 
  //              [  .   ]
  //               row j(k)            col j(k)
  //
  //     where: y_k = a * op(P2') * op(M') * e(i(k)) <: R^m
  // 
  // Of course above we only need to set the row and column elements for S(j(k),:) and S(:,j(k))
  // for the portion of the symmetric S that is being stored.
  //
  const size_type
    M_rows  = this->rows(),
    M_cols  = this->cols(),
    P1_cols = cols( P1_in.rows(), P1_in.cols(), P1_trans );
  DenseLinAlgPack::MtM_assert_sizes(
    M_rows, M_cols, trans_not(M_trans_in)
    , P1_in.rows(), P1_in.cols(), P1_trans );
  DenseLinAlgPack::MtM_assert_sizes(
    M_rows, M_cols, M_trans_in
    , P2_in.rows(), P2_in.cols(), P2_trans );
  DenseLinAlgPack::Mp_M_assert_sizes(
    S->rows(), S->cols(), no_trans
    , P1_cols, P1_cols, no_trans );
  // Rename P1 and P2 so that op(P1).rows() >= op(P2).rows()
  const bool
    reorder_P1_P2 = ( rows( P1_in.rows(), P1_in.cols(), P1_trans ) 
              < rows( P2_in.rows(), P2_in.cols(), P2_trans ) );
  const GenPermMatrixSlice
    &P1 = reorder_P1_P2 ? P2_in : P1_in,
    &P2 = reorder_P1_P2 ? P1_in : P2_in;
  const BLAS_Cpp::Transp
    M_trans  = reorder_P1_P2 ? trans_not(M_trans_in) : M_trans_in;
  // Set rows and columns of S
  const size_type
    m  = S->rows(),
    n = rows( P1.rows(), P1.cols(), P1_trans );
  DVector y_k_store(m); // ToDo: use workspace allocator!
  DVectorSlice y_k = y_k_store();
  for( GenPermMatrixSlice::const_iterator P1_itr = P1.begin(); P1_itr != P1.end(); ++P1_itr )
  {
    const size_type
      i_k = P1_trans == no_trans ? P1_itr->row_i() : P1_itr->col_j(),
      j_k = P1_trans == no_trans ? P1_itr->col_j() : P1_itr->row_i();
    // e(i(k))
    EtaVector
      e_i_k(i_k,n);
    // y_k = op(P2')*op(M')*e(i(k))
    AbstractLinAlgPack::Vp_StPtMtV( &y_k, 1.0, P2, trans_not(P2_trans), *this, trans_not(M_trans), e_i_k(), 0.0 );
    // S(j(k),:) += a*y_k'
    if( S->uplo() == BLAS_Cpp::upper )
      DenseLinAlgPack::Vp_StV( &S->gms().row(j_k)(1,j_k), a, y_k(1,j_k) );
    else
      DenseLinAlgPack::Vp_StV( &S->gms().row(j_k)(j_k,m), a, y_k(j_k,m) );
    // S(:,j(k)) += a*y_k
    if( S->uplo() == BLAS_Cpp::upper )
      DenseLinAlgPack::Vp_StV( &S->gms().col(j_k)(1,j_k), a, y_k(1,j_k) );
    else
      DenseLinAlgPack::Vp_StV( &S->gms().col(j_k)(j_k,m), a, y_k(j_k,m) );
  }
}

// Level-3 BLAS

void MatrixOpSerial::Mp_StMtM(DMatrixSlice* C, value_type a
  , BLAS_Cpp::Transp A_trans, const DMatrixSlice& B
  , BLAS_Cpp::Transp B_trans, value_type b) const
{
  DenseLinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
    , rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
  //
  // C = b*C + a*op(A)*op(B)
  // 
  // C.col(j) = b*col(j) + a*op(A)*op(B).col(j)
  //

  // ToDo: Add the ability to also perform by row if faster

  for( size_type j = 1; j <= C->cols(); ++j )
    AbstractLinAlgPack::Vp_StMtV( &C->col(j), a, *this, A_trans, DenseLinAlgPack::col( B, B_trans, j ), b );
}

void MatrixOpSerial::Mp_StMtM(DMatrixSlice* C, value_type a, const DMatrixSlice& A
  , BLAS_Cpp::Transp A_trans, BLAS_Cpp::Transp B_trans, value_type b) const
{
  DenseLinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
    , A.rows(), A.cols(), A_trans, rows(), cols(), B_trans );
  //
  // C = b*C + a*op(A)*op(B)
  // 
  // C.row(i) = b*row(i) + a*op(B)'*op(A).row(i)
  //

  // ToDo: Add the ability to also perform by column if faster

  for( size_type i = 1; i <= C->rows(); ++i )
    AbstractLinAlgPack::Vp_StMtV( &C->row(i), a, *this, BLAS_Cpp::trans_not(A_trans)
      , DenseLinAlgPack::row(A,A_trans,i) , b );
}

void MatrixOpSerial::Mp_StMtM(DMatrixSlice* C, value_type a
  , BLAS_Cpp::Transp A_trans, const MatrixOpSerial& B
  , BLAS_Cpp::Transp B_trans, value_type b) const
{
  using LinAlgOpPack::assign;
  // C = b*C + a*op(A)*op(B)
  DenseLinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
    , rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
  // Convert one of the matrices to dense, which ever one is the smallest!
  const size_type
    size_A = rows() * cols(),
    size_B = B.rows() * B.cols();
  if( size_A < size_B ) {
    DMatrix A_dense;
    assign( &A_dense, *this, BLAS_Cpp::no_trans );
    AbstractLinAlgPack::Mp_StMtM( C, a, A_dense(), A_trans, B, B_trans, b );
  }
  else {
    DMatrix B_dense;
    assign( &B_dense, B, BLAS_Cpp::no_trans );
    AbstractLinAlgPack::Mp_StMtM( C, a, *this, A_trans, B_dense(), B_trans, b );
  }
}

void MatrixOpSerial::Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceSym& sym_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // Todo: Implement!
}

void MatrixOpSerial::Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceSym& sym_rhs1
  , BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // Todo: Implement!
}

void MatrixOpSerial::Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceTri& tri_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // Todo: Implement!
}

void MatrixOpSerial::Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceTri& tri_rhs1
  , BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // Todo: Implement!
}

void MatrixOpSerial:: syrk(
  BLAS_Cpp::Transp M_trans, value_type a
  , value_type b, DMatrixSliceSym* S ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::trans_not;
  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  //
  // S = b*S + a*op(M)*op(M')
  //
  const size_type
    M_rows = this->rows(),
    M_cols = this->cols(),
    opM_rows = rows( M_rows, M_cols, M_trans ),
    opM_cols = cols( M_rows, M_cols, M_trans ),
    m = opM_rows;
  DenseLinAlgPack::Mp_MtM_assert_sizes(
    S->rows(), S->cols(), no_trans
    ,M_rows, M_cols, M_trans
    ,M_rows, M_cols, trans_not(M_trans) );
  // S = b*S
  DenseLinAlgPack::Mt_S( &DenseLinAlgPack::nonconst_tri_ele(S->gms(),S->uplo()), b );
  //
  // Form this matrix one column at a time by multiplying by e(j):
  //
  // S(:,j) += a*op(M)*(op(M')*e(j))
  //
  //    j = 1 ... opM_rows
  //
  Workspace<value_type> t1_ws(wss,opM_cols),
                           t2_ws(wss,opM_rows);
  DVectorSlice                t1(&t1_ws[0],t1_ws.size()),
                           t2(&t2_ws[0],t2_ws.size());
  for( size_type j = 1; j <= opM_rows; ++j ) {
    EtaVector e_j(j,opM_rows);
    LinAlgOpPack::V_MtV(&t1,*this,trans_not(M_trans),e_j()); // t1 = op(M')*e(j)
    LinAlgOpPack::V_MtV(&t2,*this,M_trans,t1);               // t2 = op(M)*t1
    // S(j,:) += a*t2' 
    if( S->uplo() == BLAS_Cpp::upper )
      DenseLinAlgPack::Vp_StV( &S->gms().row(j)(1,j), a, t2(1,j) );
    else
      DenseLinAlgPack::Vp_StV( &S->gms().row(j)(j,m), a, t2(j,m) );
    // S(:,j) += a*t2
    if( S->uplo() == BLAS_Cpp::upper )
      DenseLinAlgPack::Vp_StV( &S->gms().col(j)(1,j), a, t2(1,j) );
    else
      DenseLinAlgPack::Vp_StV( &S->gms().col(j)(j,m), a, t2(j,m) );
  }
}

// Overridden from MatrixOp
  
const VectorSpace& MatrixOpSerial::space_cols() const
{
  const size_type rows = this->rows();
  if(space_cols_.dim() != rows)
    space_cols_.initialize(rows);
  return space_cols_;
}
  
const VectorSpace& MatrixOpSerial::space_rows() const
{
  const size_type cols = this->cols();
  if(space_rows_.dim() != cols)
    space_rows_.initialize(cols);
  return space_rows_;
}
  
std::ostream& MatrixOpSerial::output(std::ostream& out) const {
  DMatrix tmp( 0.0, rows(), cols() );
  this->Mp_StM( &tmp(), 1.0 , BLAS_Cpp::no_trans );
  return out << tmp();
}

bool MatrixOpSerial::Mp_StM(
  MatrixOp* mwo_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs
  ) const
{
  MatrixOpGetGMSMutable
    *mwo_gms_lhs = dynamic_cast<MatrixOpGetGMSMutable*>(mwo_lhs);
  if(!mwo_gms_lhs)
    return MatrixOp::Mp_StM(mwo_lhs,alpha,trans_rhs); // boot it!
  this->Mp_StM( &MatrixDenseMutableEncap(mwo_gms_lhs)(), alpha, trans_rhs );
  return true;
}
  
bool MatrixOpSerial::Mp_StMtP(
  MatrixOp* mwo_lhs, value_type alpha
  , BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  ) const
{
  MatrixOpGetGMSMutable
    *mwo_gms_lhs = dynamic_cast<MatrixOpGetGMSMutable*>(mwo_lhs);
  if(!mwo_gms_lhs)
    return MatrixOp::Mp_StMtP(mwo_lhs,alpha,M_trans,P_rhs,P_rhs_trans); // boot it!
  this->Mp_StMtP(&MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,M_trans,P_rhs,P_rhs_trans);
  return true;
}
  
bool MatrixOpSerial::Mp_StPtM(
  MatrixOp* mwo_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  , BLAS_Cpp::Transp M_trans
  ) const
{
  MatrixOpGetGMSMutable
    *mwo_gms_lhs = dynamic_cast<MatrixOpGetGMSMutable*>(mwo_lhs);
  if(!mwo_gms_lhs)
    return MatrixOp::Mp_StPtM(mwo_lhs,alpha,P_rhs,P_rhs_trans,M_trans); // boot it!
  this->Mp_StPtM(&MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,P_rhs,P_rhs_trans,M_trans);
  return true;
}
  
bool MatrixOpSerial::Mp_StPtMtP(
  MatrixOp* mwo_lhs, value_type alpha
  ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  ,BLAS_Cpp::Transp M_trans
  ,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  ) const
{
  MatrixOpGetGMSMutable
    *mwo_gms_lhs = dynamic_cast<MatrixOpGetGMSMutable*>(mwo_lhs);
  if(!mwo_gms_lhs)
    return MatrixOp::Mp_StPtMtP(mwo_lhs,alpha,P_rhs1,P_rhs1_trans,M_trans,P_rhs2,P_rhs2_trans); // boot it!
  this->Mp_StPtMtP(&MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,P_rhs1,P_rhs1_trans,M_trans,P_rhs2,P_rhs2_trans);
  return true;
}
  
void MatrixOpSerial::Vp_StMtV(
  VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const Vector& v_rhs2, value_type beta) const
{
  VectorDenseMutableEncap       vs_lhs(*v_lhs);
  VectorDenseEncap              vs_rhs2(v_rhs2);
  this->Vp_StMtV( &vs_lhs(), alpha, trans_rhs1, vs_rhs2(), beta );	
}

void MatrixOpSerial::Vp_StMtV(
  VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2, value_type beta) const
{
  VectorDenseMutableEncap       vs_lhs(*v_lhs);
  this->Vp_StMtV( &vs_lhs(), alpha, trans_rhs1, sv_rhs2, beta );	
}
  
void MatrixOpSerial::Vp_StPtMtV(
  VectorMutable* v_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const Vector& v_rhs3, value_type beta) const
{
  VectorDenseMutableEncap       vs_lhs(*v_lhs);
  VectorDenseEncap              vs_rhs3(v_rhs3);
  this->Vp_StPtMtV( &vs_lhs(), alpha, P_rhs1, P_rhs1_trans, M_rhs2_trans, vs_rhs3(), beta );	
}
  
void MatrixOpSerial::Vp_StPtMtV(
  VectorMutable* v_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , BLAS_Cpp::Transp M_rhs2_trans
  , const SpVectorSlice& sv_rhs3, value_type beta) const
{
  VectorDenseMutableEncap       vs_lhs(*v_lhs);
  this->Vp_StPtMtV( &vs_lhs(), alpha, P_rhs1, P_rhs1_trans, M_rhs2_trans, sv_rhs3, beta );	
}
  
value_type MatrixOpSerial::transVtMtV(
  const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const Vector& v_rhs3) const
{
  VectorDenseEncap              vs_rhs1(v_rhs1);
  VectorDenseEncap              vs_rhs3(v_rhs3);
  return this->transVtMtV(vs_rhs1(),trans_rhs2,vs_rhs3());
}
  
void MatrixOpSerial::syr2k(
  BLAS_Cpp::Transp M_trans, value_type alpha
  , const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
  , const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
  , value_type beta, MatrixSymOp* symwo_lhs ) const
{
  MatrixSymOpGetGMSSymMutable
    *symwo_gms_lhs = dynamic_cast<MatrixSymOpGetGMSSymMutable*>(symwo_lhs);
  if(!symwo_gms_lhs) {
    MatrixOp::syr2k(M_trans,alpha,P1,P1_trans,P2,P2_trans,beta,symwo_lhs); // Boot it
    return;
  }
  this->syr2k(
    M_trans,alpha,P1,P1_trans,P2,P2_trans,beta
    ,&MatrixDenseSymMutableEncap(symwo_gms_lhs)()
     );
}

bool MatrixOpSerial::Mp_StMtM(
  MatrixOp* mwo_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOp& mwo_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta ) const
{
  MatrixOpGetGMSMutable
    *mwo_gms_lhs = dynamic_cast<MatrixOpGetGMSMutable*>(mwo_lhs);
  if(mwo_gms_lhs) {
    // Try to match the rhs arguments to known serial interfaces
    if(const MatrixOpGetGMS* mwo_gms_rhs2 = dynamic_cast<const MatrixOpGetGMS*>(&mwo_rhs2)) {
      this->Mp_StMtM(
        &MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,trans_rhs1
        ,MatrixDenseEncap(*mwo_gms_rhs2)(),trans_rhs2,beta );
      return true;
    }
    if(const MatrixSymOpGetGMSSym* mwo_sym_gms_rhs2 = dynamic_cast<const MatrixSymOpGetGMSSym*>(&mwo_rhs2)) {
      this->Mp_StMtM(
        &MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,trans_rhs1
        ,MatrixDenseEncap(*mwo_sym_gms_rhs2)(),trans_rhs2,beta );
      return true;
    }
    if(const MatrixOpGetGMSTri* mwo_tri_gms_rhs2 = dynamic_cast<const MatrixOpGetGMSTri*>(&mwo_rhs2)) {
      this->Mp_StMtM(
        &MatrixDenseMutableEncap(mwo_gms_lhs)(),alpha,trans_rhs1
        ,MatrixDenseEncap(*mwo_tri_gms_rhs2)(),trans_rhs2,beta );
      return true;
    }
    // If we get here, the matrix arguments did not match up so we have to give up (I think?)
  }
  // Let the default implementation try to find matrix arguments that can handle this!
  return MatrixOp::Mp_StMtM(mwo_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2,beta); // Boot it!
}

bool MatrixOpSerial::syrk(
  BLAS_Cpp::Transp M_trans, value_type alpha
  , value_type beta, MatrixSymOp* symwo_lhs ) const
{
  MatrixSymOpGetGMSSymMutable
    *symwo_gms_lhs = dynamic_cast<MatrixSymOpGetGMSSymMutable*>(symwo_lhs);
  if(!symwo_gms_lhs) {
    return MatrixOp::syrk(M_trans,alpha,beta,symwo_lhs); // Boot it
  }
  this->syrk(M_trans,alpha,beta,&MatrixDenseSymMutableEncap(symwo_gms_lhs)());
  return true;
}

}	// end namespace AbstractLinAlgPack
