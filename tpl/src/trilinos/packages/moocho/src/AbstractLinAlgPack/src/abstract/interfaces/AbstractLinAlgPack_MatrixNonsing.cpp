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

#include "AbstractLinAlgPack_MatrixNonsing.hpp"
#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorView.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Clone

MatrixNonsing::mat_mns_mut_ptr_t
MatrixNonsing::clone_mns()
{
  return Teuchos::null;
}

MatrixNonsing::mat_mns_ptr_t
MatrixNonsing::clone_mns() const
{
  return const_cast<MatrixNonsing*>(this)->clone_mns(); // Implicit conversion to const
}

// Level-2 BLAS

void MatrixNonsing::V_InvMtV(
  VectorMutable* y, BLAS_Cpp::Transp M_trans, const SpVectorSlice& sx
  ) const
{
  if( sx.nz() ) {
    VectorSpace::vec_mut_ptr_t
      x = (M_trans == BLAS_Cpp::no_trans
            ? this->space_cols()
            : this->space_rows()
        ).create_member();
    x->set_sub_vector(sub_vec_view(sx));
    this->V_InvMtV(y,M_trans,*x);
  }
  else {
    *y = 0.0;
  }
}

value_type MatrixNonsing::transVtInvMtV(
  const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2, const Vector& v_rhs3
  ) const
{
  VectorSpace::vec_mut_ptr_t
    v = (trans_rhs2 == BLAS_Cpp::no_trans
          ? this->space_rows()
          : this->space_cols()
      ).create_member();
  this->V_InvMtV( v.get(), trans_rhs2, v_rhs3 );
  return dot(v_rhs1,*v);
}

value_type MatrixNonsing::transVtInvMtV(
  const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3
  ) const
{
  VectorSpace::vec_mut_ptr_t
    v = (trans_rhs2 == BLAS_Cpp::no_trans
          ? this->space_rows()
          : this->space_cols()
      ).create_member();
  this->V_InvMtV( v.get(), trans_rhs2, sv_rhs3 );
  return dot(sv_rhs1,*v);
}

// Level-3 BLAS

void MatrixNonsing::M_StInvMtM(
  MatrixOp* C_lhs, value_type alpha
  ,BLAS_Cpp::Transp M_trans
  ,const MatrixOp& B, BLAS_Cpp::Transp B_trans
  ) const
{
  //
  // C = a * inv(op(M)) * op(B)
  //
  using Teuchos::dyn_cast;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    C_lhs == NULL, std::invalid_argument
    ,"MatrixNonsing::M_StInvMtM(...) : Error!" );
  
#endif
  const size_type
    C_rows = C_lhs->rows(),
    C_cols = C_lhs->cols();
  const size_type
    op_B_cols = BLAS_Cpp::cols( B.rows(), B.cols(), B_trans );
#ifdef TEUCHOS_DEBUG
  // We can't check vector spaces since *this may not support MatrixOp
  // However, we could dynamic cast to see if MatrixOp is supported and then
  // be able to use Mp_MtM_assert_compatibility() but this is okay for now.
  const size_type
    M_rows    = this->rows(),
    M_cols    = this->cols(),
    op_B_rows = BLAS_Cpp::rows( B.rows(), B.cols(), B_trans );
  TEUCHOS_TEST_FOR_EXCEPTION(
    C_rows != M_rows || M_rows != M_cols || M_cols != op_B_rows || C_cols != op_B_cols
    , std::invalid_argument
    ,"MatrixNonsing::M_StInvMtM(...) : Error!" );
#endif
  //
  // Compute C = a * inv(op(M)) * op(B) one column at a time:
  //
  // C(:,j) = inv(op(M)) * a * op(B) * e(j)    , for j = 1...C.cols()
  //                       \______________/    
  //                              t_j
  //
  MultiVectorMutable  &C = dyn_cast<MultiVectorMutable>(*C_lhs);
  VectorSpace::vec_mut_ptr_t
    t_j = ( B_trans == no_trans ? B.space_cols() : B.space_rows() ).create_member();
  for( size_type j = 1; j <= C_cols; ++j ) {
    // t_j = alpha * op(B) * e_j
    EtaVector e_j( j, op_B_cols );
    LinAlgOpPack::V_StMtV( t_j.get(), alpha, B, B_trans, e_j() );
    // C(:,j) = inv(op(M)) * t_j
    AbstractLinAlgPack::V_InvMtV( C.col(j).get(), *this, M_trans, *t_j );
  }
}

void MatrixNonsing::M_StMtInvM(
  MatrixOp* g_lhs, value_type alpha
  ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,BLAS_Cpp::Transp trans_rhs2
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
}

}	// end namespace AbstractLinAlgPack
