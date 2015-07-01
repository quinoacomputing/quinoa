#if 0

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

#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "AbstractLinAlgPack_COOMatrixPartitionViewSubclass.hpp"
#include "AbstractLinAlgPack_SparseVectorSliceOp.hpp"
#include "AbstractLinAlgPack_SparseElement.hpp"
#include "AbstractLinAlgPack_COOMPartitionOp.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"

namespace LinAlgOpPack {

using AbstractLinAlgPack::Vp_StV;
using AbstractLinAlgPack::Vp_StMtV;
using AbstractLinAlgPack::Mp_StM;
using AbstractLinAlgPack::Mp_StMtM;

}	// end namespace LinAlgOpPack

namespace AbstractLinAlgPack {

size_type COOMatrixPartitionViewSubclass::rows() const {
  return trans_ == BLAS_Cpp::no_trans ? m().rows() : m().cols();
}

size_type COOMatrixPartitionViewSubclass::cols() const {
  return trans_ == BLAS_Cpp::no_trans ? m().cols() : m().rows();
}

MatrixOp& COOMatrixPartitionViewSubclass::operator=(const MatrixOp& m) {
  if(&m == this) return *this;	// assignment to self
  const COOMatrixPartitionViewSubclass *p_m = dynamic_cast<const COOMatrixPartitionViewSubclass*>(&m);
  if(p_m) {
    throw std::invalid_argument("COOMatrixPartitionViewSubclass::operator=(const MatrixOp& m)"
      " :  There is not an assignment operator defined for COOMatrixWithPartitionedView::partition_type"
      ".   Only assignment to self can be handeled" );
  }
  else {
    throw std::invalid_argument("COOMatrixPartitionViewSubclass::operator=(const MatrixOp& m)"
      " : The concrete type of m is not a subclass of COOMatrixPartitionViewSubclass as expected" );
  }
  return *this;
}

// Level-1 BLAS

void COOMatrixPartitionViewSubclass::Mp_StM(DMatrixSlice* gms_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs) const
{
  AbstractLinAlgPack::Mp_StM(gms_lhs,alpha,m(),op(trans_rhs));
}

// Level-2 BLAS

void COOMatrixPartitionViewSubclass::Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2, value_type beta) const
{
  AbstractLinAlgPack::Vp_StMtV(vs_lhs, alpha, m(), op(trans_rhs1), vs_rhs2, beta);
}

void COOMatrixPartitionViewSubclass::Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2, value_type beta) const
{
  DVector v_rhs2;
  LinAlgOpPack::assign(&v_rhs2,sv_rhs2);
  AbstractLinAlgPack::Vp_StMtV(vs_lhs, alpha, m(), op(trans_rhs1), v_rhs2(), beta);
}

value_type COOMatrixPartitionViewSubclass::transVtMtV(const DVectorSlice& vs_rhs1
  , BLAS_Cpp::Transp trans_rhs2, const DVectorSlice& vs_rhs3) const
{
  DVector tmp;
  LinAlgOpPack::V_MtV(&tmp,m(),op(trans_rhs2),vs_rhs3);
  return DenseLinAlgPack::dot(vs_rhs1,tmp());
}

value_type COOMatrixPartitionViewSubclass::transVtMtV(const SpVectorSlice& sv_rhs1
  , BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const
{
  DVector v_rhs3;
  LinAlgOpPack::assign(&v_rhs3,sv_rhs3);
  DVector tmp;
  LinAlgOpPack::V_MtV(&tmp,m(),op(trans_rhs2),v_rhs3());
  return dot(sv_rhs1,tmp());
}

// Level-3 BLAS

void COOMatrixPartitionViewSubclass::Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha
  , BLAS_Cpp::Transp trans_rhs1, const DMatrixSlice& gms_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
  AbstractLinAlgPack::Mp_StMtM(gms_lhs, alpha, m(), op(trans_rhs1), gms_rhs2, trans_rhs2, beta);
}

void COOMatrixPartitionViewSubclass::Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSlice& gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const
{
  AbstractLinAlgPack::Mp_StMtM(gms_lhs, alpha, gms_rhs1, trans_rhs1, m(), op(trans_rhs2), beta);
}

}	// end namespace AbstractLinAlgPack

#endif // 0
