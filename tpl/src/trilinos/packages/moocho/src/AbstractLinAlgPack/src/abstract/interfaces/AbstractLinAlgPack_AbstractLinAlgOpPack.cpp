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

#include "AbstractLinAlgPack_LinAlgOpPack.hpp"

// Level 1 BLAS for Matrices

// M_lhs = op(M_rhs).

void LinAlgOpPack::assign(MatrixOp* M_lhs, const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
  Mp_M_assert_compatibility( M_lhs, BLAS_Cpp::no_trans, M_rhs, trans_rhs );
  M_lhs->zero_out();
  Mp_StM(M_lhs,1.0,M_rhs,trans_rhs);
}

// M_lhs = alpha * op(M_rhs).

void LinAlgOpPack::M_StM(MatrixOp* M_lhs, value_type alpha, const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
  Mp_M_assert_compatibility( M_lhs, BLAS_Cpp::no_trans, M_rhs, trans_rhs );
  M_lhs->zero_out();
  Mp_StM(M_lhs,alpha,M_rhs,trans_rhs);
}

// M_lhs = op(M_rhs1) + op(M_rhs2).

void LinAlgOpPack::M_MpM(MatrixOp* M_lhs, const MatrixOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_M_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1);
  MopM_assert_compatibility(M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  M_lhs->zero_out();
  Mp_M(M_lhs,M_rhs1,trans_rhs1);
  Mp_M(M_lhs,M_rhs2,trans_rhs2);
}

// M_lhs = op(M_rhs1) - op(M_rhs2).

void LinAlgOpPack::M_MmM(MatrixOp* M_lhs, const MatrixOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_M_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1);
  MopM_assert_compatibility(M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  M_lhs->zero_out();
  Mp_M(M_lhs,M_rhs1,trans_rhs1);
  Mp_StM(M_lhs,-1.0,M_rhs2,trans_rhs2);
}

// M_lhs = alpha * op(M_rhs1) + op(m_rhs2).

void LinAlgOpPack::M_StMpM(MatrixOp* M_lhs, value_type alpha, const MatrixOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_M_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1);
  MopM_assert_compatibility(M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  assign(M_lhs,M_rhs2,trans_rhs2);
  Mp_StM(M_lhs,alpha,M_rhs1,trans_rhs1);
}

// Level 3 BLAS

// M_lhs = alpha * op(M_rhs1) * op(M_rhs2).

void LinAlgOpPack::M_StMtM(MatrixOp* M_lhs, value_type alpha, const MatrixOp& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_MtM_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  Mp_StMtM(M_lhs,alpha,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2,0.0);
}

// M_lhs = op(M_rhs1) * op(M_rhs2).

void LinAlgOpPack::M_MtM(MatrixOp* M_lhs, const MatrixOp& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
  Mp_MtM_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
  Mp_StMtM(M_lhs,1.0,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2,0.0);
}
