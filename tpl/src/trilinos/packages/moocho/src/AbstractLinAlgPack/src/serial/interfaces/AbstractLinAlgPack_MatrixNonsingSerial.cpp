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

// ToDo: 3/6/00: Provide default implementations for these
// operations.

#include <assert.h>

#include "AbstractLinAlgPack_MatrixNonsingSerial.hpp"
#include "AbstractLinAlgPack_MatrixOpSerial.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_MatrixOpGetGMSMutable.hpp"
#include "AbstractLinAlgPack_MatrixOpGetGMSTri.hpp"
#include "AbstractLinAlgPack_MatrixSymOpGetGMSSymMutable.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::Mp_StM;
  using AbstractLinAlgPack::Vp_StMtV;
}

namespace AbstractLinAlgPack {

//  Level-2 BLAS

void MatrixNonsingSerial::V_InvMtV(
  DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1,const DVectorSlice& vs_rhs2
  ) const
{
  const size_type n = rows();
  DenseLinAlgPack::MtV_assert_sizes( n, n, trans_rhs1, vs_rhs2.dim() );
  v_lhs->resize(n);
  this->V_InvMtV( &(*v_lhs)(), trans_rhs1, vs_rhs2 );
}

void MatrixNonsingSerial::V_InvMtV(
  DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2
  ) const
{
  const size_type n = rows();
  DenseLinAlgPack::MtV_assert_sizes( n, n, trans_rhs1, sv_rhs2.dim() );
  v_lhs->resize(n);
  DVector v_rhs2;
  LinAlgOpPack::assign( &v_rhs2, sv_rhs2 );
  this->V_InvMtV( &(*v_lhs)(), trans_rhs1, v_rhs2() );
}

void MatrixNonsingSerial::V_InvMtV(
  DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2
  ) const
{
  const size_type n = rows();
  DenseLinAlgPack::Vp_MtV_assert_sizes( vs_lhs->dim(), n, n, trans_rhs1, sv_rhs2.dim() );
  DVector v_rhs2;
  LinAlgOpPack::assign( &v_rhs2, sv_rhs2 );
  this->V_InvMtV( vs_lhs, trans_rhs1, v_rhs2() );
}

value_type MatrixNonsingSerial::transVtInvMtV(
  const DVectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2, const DVectorSlice& vs_rhs3
  ) const
{
  const size_type n = rows();
  DenseLinAlgPack::Vp_MtV_assert_sizes( vs_rhs1.dim(), n, n, trans_rhs2, vs_rhs3.dim() );
  DVector tmp;
  this->V_InvMtV( &tmp, trans_rhs2, vs_rhs3 );
  return DenseLinAlgPack::dot( vs_rhs1, tmp() );
}

value_type MatrixNonsingSerial::transVtInvMtV(
  const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3
  ) const
{
  const size_type n = rows();
  DenseLinAlgPack::Vp_MtV_assert_sizes( sv_rhs1.dim(), n, n, trans_rhs2, sv_rhs3.dim() );
  DVector tmp;
  this->V_InvMtV( &tmp, trans_rhs2, sv_rhs3 );
  return AbstractLinAlgPack::dot( sv_rhs1, tmp() );
}

// Level-3 BLAS

void MatrixNonsingSerial::M_StInvMtM(
  DMatrix* C, value_type a
  ,BLAS_Cpp::Transp A_trans
  ,const DMatrixSlice& B, BLAS_Cpp::Transp B_trans
  ) const
{
  DenseLinAlgPack::MtM_assert_sizes( rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
  C->resize(
      BLAS_Cpp::rows( rows(), cols(), A_trans )
    , BLAS_Cpp::cols( B.rows(), B.cols(), B_trans )
    );
  this->M_StInvMtM( &(*C)(), a, A_trans, B, B_trans );
}

void MatrixNonsingSerial::M_StInvMtM(
  DMatrixSlice* C, value_type a
  ,BLAS_Cpp::Transp A_trans
  ,const DMatrixSlice& B, BLAS_Cpp::Transp B_trans
  ) const
{
  DenseLinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
    , rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
  //
  // C = a * inv(op(A)) * op(B)
  //
  // C.col(j) = a * inv(op(A)) * op(B).col(j)
  //

  for( size_type j = 1; j <= C->cols(); ++j )
    AbstractLinAlgPack::V_InvMtV( &C->col(j), *this, A_trans
      , DenseLinAlgPack::col( B, B_trans, j ) );
  if( a != 1.0 )
    LinAlgOpPack::Mt_S( C, a );
}

void MatrixNonsingSerial::M_StMtInvM(
  DMatrix* gm_lhs, value_type alpha
  ,const DMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);	// ToDo: Implement this!
}

void MatrixNonsingSerial::M_StMtInvM(
  DMatrixSlice* gms_lhs, value_type alpha
  ,const DMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);	// ToDo: Implement this!
}

void MatrixNonsingSerial::M_StInvMtM(
  DMatrix* C, value_type a
  ,BLAS_Cpp::Transp A_trans
  ,const MatrixOpSerial& B, BLAS_Cpp::Transp B_trans
  ) const
{
  DenseLinAlgPack::MtM_assert_sizes( rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
  C->resize(
      BLAS_Cpp::rows( rows(), cols(), A_trans )
    , BLAS_Cpp::cols( B.rows(), B.cols(), B_trans )
    );
  AbstractLinAlgPack::M_StInvMtM( &(*C)(), a, *this, A_trans, B, B_trans );
}

void MatrixNonsingSerial::M_StInvMtM(
  DMatrixSlice* C, value_type a
  ,BLAS_Cpp::Transp A_trans
  ,const MatrixOpSerial& B, BLAS_Cpp::Transp B_trans
  ) const
{
  using LinAlgOpPack::assign;
  DenseLinAlgPack::Mp_MtM_assert_sizes( C->rows(), C->cols(), BLAS_Cpp::no_trans
    , rows(), cols(), A_trans, B.rows(), B.cols(), B_trans );
  DMatrix B_dense;
  assign( &B_dense, B, BLAS_Cpp::no_trans );
  AbstractLinAlgPack::M_StInvMtM( C, a, *this, A_trans, B_dense(), B_trans );
}

void MatrixNonsingSerial::M_StMtInvM(
  DMatrix* gm_lhs, value_type alpha
  ,const MatrixOpSerial& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);	// ToDo: Implement this!
}

void MatrixNonsingSerial::M_StMtInvM(
  DMatrixSlice* gms_lhs, value_type alpha
  ,const MatrixOpSerial& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);	// ToDo: Implement this!
}

// Overridden from MatrixNonsing

void MatrixNonsingSerial::V_InvMtV(
  VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
  ,const Vector& v_rhs2) const
{
  VectorDenseMutableEncap       vs_lhs(*v_lhs);
  VectorDenseEncap              vs_rhs2(v_rhs2);
  this->V_InvMtV( &vs_lhs(), trans_rhs1, vs_rhs2() );	
}

void MatrixNonsingSerial::V_InvMtV(
  VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
  ,const SpVectorSlice& sv_rhs2) const
{
  this->V_InvMtV( &VectorDenseMutableEncap(*v_lhs)(), trans_rhs1, sv_rhs2 );
}

value_type MatrixNonsingSerial::transVtInvMtV(
  const Vector& v_rhs1
  ,BLAS_Cpp::Transp trans_rhs2, const Vector& v_rhs3) const
{
  VectorDenseEncap              vs_rhs1(v_rhs1);
  VectorDenseEncap              vs_rhs3(v_rhs3);
  return this->transVtInvMtV(vs_rhs1(),trans_rhs2,vs_rhs3());
}

void MatrixNonsingSerial::M_StInvMtM(
  MatrixOp* m_lhs, value_type alpha
  ,BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
  ) const
{
  using Teuchos::dyn_cast;
  MatrixDenseMutableEncap
    gms_lhs(m_lhs);      // Warning!  This may throw an exception!
  if(const MatrixOpGetGMS* mwo_gms_rhs2 = dynamic_cast<const MatrixOpGetGMS*>(&mwo_rhs2)) {
    this->M_StInvMtM(&gms_lhs(),alpha,trans_rhs1,MatrixDenseEncap(*mwo_gms_rhs2)(),trans_rhs2);
    return;
  }
  this->M_StInvMtM(&gms_lhs(),alpha,trans_rhs1,dyn_cast<const MatrixOpSerial>(mwo_rhs2),trans_rhs2);
}

void MatrixNonsingSerial::M_StMtInvM(
  MatrixOp* m_lhs, value_type alpha
  ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2
  ) const
{
  using Teuchos::dyn_cast;
  MatrixDenseMutableEncap
    gms_lhs(m_lhs);      // Warning!  This may throw an exception!
  if(const MatrixOpGetGMS* mwo_gms_rhs1 = dynamic_cast<const MatrixOpGetGMS*>(&mwo_rhs1)) {
    this->M_StMtInvM(&gms_lhs(),alpha,MatrixDenseEncap(*mwo_gms_rhs1)(),trans_rhs1,trans_rhs2);
    return;
  }
  this->M_StMtInvM(&gms_lhs(),alpha,dyn_cast<const MatrixOpSerial>(mwo_rhs1),trans_rhs1,trans_rhs2);
}

}	// end namespace AbstractLinAlgPack
