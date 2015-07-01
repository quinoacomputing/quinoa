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
//

#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "AbstractLinAlgPack_VectorMutableDense.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_MatrixOpGetGMS.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"


void LinAlgOpPack::assign(
  DMatrixSlice* gms_lhs, const MatrixOp& M_rhs,
  BLAS_Cpp::Transp trans_rhs
  )
{
  Mp_M_assert_sizes(
    gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans,
    M_rhs.rows(), M_rhs.cols(), trans_rhs );
  (*gms_lhs) = 0.0;
  Mp_StM(gms_lhs,1.0,M_rhs,trans_rhs);
}


void LinAlgOpPack::Mp_StM(
  DMatrixSlice* C, value_type a
  ,const MatrixOp& B, BLAS_Cpp::Transp B_trans
  )
{
  using AbstractLinAlgPack::VectorSpace;
  using AbstractLinAlgPack::VectorDenseEncap;
  using AbstractLinAlgPack::MatrixOpGetGMS;
  using AbstractLinAlgPack::MatrixDenseEncap;
  const MatrixOpGetGMS
    *B_get_gms = dynamic_cast<const MatrixOpGetGMS*>(&B);
  if(B_get_gms) {
    DenseLinAlgPack::Mp_StM( C, a, MatrixDenseEncap(*B_get_gms)(), B_trans );		
  }
  else {
    const size_type num_cols = C->cols();
    VectorSpace::multi_vec_mut_ptr_t
      B_mv = ( B_trans == BLAS_Cpp::no_trans 
          ? B.space_cols() : B.space_rows()
          ).create_members(num_cols);
    assign(B_mv.get(),B,B_trans);
    for( size_type j = 1; j <= num_cols; ++j ) {
      DenseLinAlgPack::Vp_StV(&C->col(j),a,VectorDenseEncap(*B_mv->col(j))());
    }
  }
}

void LinAlgOpPack::Vp_StMtV(
  DVectorSlice* y, value_type a, const MatrixOp& M
  ,BLAS_Cpp::Transp M_trans, const DVectorSlice& x, value_type b
  )
{
  using BLAS_Cpp::no_trans;
  using AbstractLinAlgPack::VectorDenseMutableEncap;
  VectorSpace::vec_mut_ptr_t
    ay = ( M_trans == no_trans ? M.space_cols() : M.space_rows() ).create_member(),
    ax = ( M_trans == no_trans ? M.space_rows() : M.space_cols() ).create_member();
  (VectorDenseMutableEncap(*ay))() = *y;
  (VectorDenseMutableEncap(*ax))() = x;
  AbstractLinAlgPack::Vp_StMtV( ay.get(), a, M, M_trans, *ax, b );
  *y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::Vp_StMtV(
  DVectorSlice* y, value_type a, const MatrixOp& M
  ,BLAS_Cpp::Transp M_trans, const SpVectorSlice& x, value_type b
  )
{
  using BLAS_Cpp::no_trans;
  using AbstractLinAlgPack::VectorDenseMutableEncap;
  VectorSpace::vec_mut_ptr_t
    ay = ( M_trans == no_trans ? M.space_cols() : M.space_rows() ).create_member();
  (VectorDenseMutableEncap(*ay))() = *y;
  AbstractLinAlgPack::Vp_StMtV( ay.get(), a, M, M_trans, x, b );
  *y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::V_InvMtV(
  DVectorSlice* y, const MatrixOpNonsing& M
  ,BLAS_Cpp::Transp M_trans, const DVectorSlice& x
  )
{
  using BLAS_Cpp::trans;
  using AbstractLinAlgPack::VectorDenseMutableEncap;
  VectorSpace::vec_mut_ptr_t
    ay = ( M_trans == trans ? M.space_cols() : M.space_rows() ).create_member(),
    ax = ( M_trans == trans ? M.space_rows() : M.space_cols() ).create_member();
  (VectorDenseMutableEncap(*ay))() = *y;
  (VectorDenseMutableEncap(*ax))() = x;
  AbstractLinAlgPack::V_InvMtV( ay.get(), M, M_trans, *ax );
  *y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::V_InvMtV(
  DVector* y, const MatrixOpNonsing& M
  ,BLAS_Cpp::Transp M_trans, const DVectorSlice& x
  )
{
  using BLAS_Cpp::trans;
  using AbstractLinAlgPack::VectorDenseMutableEncap;
  VectorSpace::vec_mut_ptr_t
    ay = ( M_trans == trans ? M.space_cols() : M.space_rows() ).create_member(),
    ax = ( M_trans == trans ? M.space_rows() : M.space_cols() ).create_member();
  (VectorDenseMutableEncap(*ax))() = x;
  AbstractLinAlgPack::V_InvMtV( ay.get(), M, M_trans, *ax );
  *y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::V_InvMtV(
  DVectorSlice* y, const MatrixOpNonsing& M
  ,BLAS_Cpp::Transp M_trans, const SpVectorSlice& x
  )
{
  using BLAS_Cpp::trans;
  using AbstractLinAlgPack::VectorDenseMutableEncap;
  VectorSpace::vec_mut_ptr_t
    ay = ( M_trans == trans ? M.space_cols() : M.space_rows() ).create_member();
  AbstractLinAlgPack::V_InvMtV( ay.get(), M, M_trans, x );
  *y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::V_InvMtV(
  DVector* y, const MatrixOpNonsing& M
  ,BLAS_Cpp::Transp M_trans, const SpVectorSlice& x
  )
{
  y->resize(M.rows());
  LinAlgOpPack::V_InvMtV( &(*y)(), M, M_trans, x );
}

// These methods below are a real problem to implement in general.
//
// If the column space of op(M) could not return the above vector space
// for op(M).space_cols().space(P,P_trans) then we will try, as a last
// resort, to a dense serial vector and see what happens.

void LinAlgOpPack::Vp_StPtMtV(
  DVectorSlice* y, value_type a
  ,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  ,const MatrixOp& M, BLAS_Cpp::Transp M_trans
  ,const DVectorSlice& x, value_type b
  )
{
  namespace mmp = MemMngPack;
  using BLAS_Cpp::no_trans;
  using AbstractLinAlgPack::VectorMutableDense;
  using AbstractLinAlgPack::VectorDenseMutableEncap;
  using AbstractLinAlgPack::Vp_StPtMtV;
  VectorSpace::space_ptr_t
    ay_space = ( M_trans == no_trans ? M.space_cols() : M.space_rows() ).space(P,P_trans);
  VectorSpace::vec_mut_ptr_t
    ay =  ( ay_space.get()
        ? ay_space->create_member()
        : Teuchos::rcp_implicit_cast<VectorMutable>(
          Teuchos::rcp(new VectorMutableDense(BLAS_Cpp::rows(P.rows(),P.cols(),P_trans)))
          ) ),
    ax = ( M_trans == no_trans ? M.space_rows() : M.space_cols() ).create_member();
  (VectorDenseMutableEncap(*ay))() = *y;
  (VectorDenseMutableEncap(*ax))() = x;
  Vp_StPtMtV( ay.get(), a, P, P_trans, M, M_trans, *ax, b );
  *y = VectorDenseMutableEncap(*ay)();
}

void LinAlgOpPack::Vp_StPtMtV(
  DVectorSlice* y, value_type a
  ,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  ,const MatrixOp& M, BLAS_Cpp::Transp M_trans
  ,const SpVectorSlice& x, value_type b
  )
{
  using BLAS_Cpp::no_trans;
  using AbstractLinAlgPack::VectorMutableDense;
  using AbstractLinAlgPack::VectorDenseMutableEncap;
  using AbstractLinAlgPack::Vp_StPtMtV;
  VectorSpace::vec_mut_ptr_t
    ay = ( M_trans == no_trans ? M.space_cols() : M.space_rows() ).space(P,P_trans)->create_member();
  (VectorDenseMutableEncap(*ay))() = *y;
  Vp_StPtMtV( ay.get(), a, P, P_trans, M, M_trans, x, b );
  *y = VectorDenseMutableEncap(*ay)();
}
