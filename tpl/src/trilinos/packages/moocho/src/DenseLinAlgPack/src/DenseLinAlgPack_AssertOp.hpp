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

#ifndef LIN_ALG_PACK_ASSERT_OP_H
#define LIN_ALG_PACK_ASSERT_OP_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

#ifdef LINALGPACK_CHECK_RHS_SIZES


/* * @name Assertion functions for linear algebra operations.
  *
  * These functions check the sizes of the linear algebra
  * expressions and throw a std::length_error if
  * the sizes do not match.  These functions
  * only perform there operations if #LINALGPACK_CHECK_RHS_SIZES#
  * is defined.
  */
// @{

/* * @name Level 1 BLAS
  */
// @{

/// v_lhs += op v_rhs
void Vp_V_assert_sizes(size_type v_lhs_size, size_type v_rhs_size);

/// v_rhs1 op v_rhs2
void VopV_assert_sizes(size_type v_rhs1_size, size_type v_rhs2_size);

/// op(m_lhs) += op op(m_rhs)
void Mp_M_assert_sizes(size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
  , size_type m_rhs_rows, size_type m_rhs_cols, BLAS_Cpp::Transp trans_rhs);

/// v_rhs1 op v_rhs2
void MopM_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2);

//		end Level 1 BLAS
// @}

/* *  @name Level 2 BLAS
  */
// @{

/// op(m_rhs1) * v_rhs2
void MtV_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols
  , BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size);

/// v_lhs += op(m_rhs1) * v_rhs2
void Vp_MtV_assert_sizes(size_type v_lhs_size, size_type m_rhs1_rows
  , size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size);

//		end Level 2 BLAS
// @}

/* *  @name Level 3 BLAS
  */
// @{

/// op(m_lhs) += op(m_rhs1)
void MtM_assert_sizes(
    size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2);

/// op(m_lhs) += op(m_rhs1) * op(m_rhs2)
void Mp_MtM_assert_sizes(
    size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
  , size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2);

//		end Level 3 BLAS
// @}


// @}

#else

// inline definitions that do nothing

inline
void Vp_V_assert_sizes(size_type v_lhs_size, size_type v_rhs_size)
{}

inline
void VopV_assert_sizes(size_type v_rhs1_size, size_type v_rhs2_size)
{}

inline
void Mp_M_assert_sizes(size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
  , size_type m_rhs_rows, size_type m_rhs_cols, BLAS_Cpp::Transp trans_rhs)
{}

inline
void MopM_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{}

inline
void MtV_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols
  , BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{}

inline
void Vp_MtV_assert_sizes(size_type v_lhs_size, size_type m_rhs1_rows
  , size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{}

inline
void MtM_assert_sizes(
    size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{}

inline
void Mp_MtM_assert_sizes(
    size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
  , size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{}

#endif

} // end namespace DenseLinAlgPack

#endif	// LIN_ALG_PACK_ASSERT_OP_H
