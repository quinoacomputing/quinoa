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

#include <stdexcept>
#include <string>

#include "DenseLinAlgPack_AssertOp.hpp"
#include "Teuchos_Assert.hpp"

#ifdef LINALGPACK_CHECK_RHS_SIZES

void DenseLinAlgPack::Vp_V_assert_sizes(size_type v_lhs_size, size_type v_rhs_size)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    v_lhs_size != v_rhs_size, std::length_error
    ,"Vp_V_assert_sizes(...) : The sizes of v_lhs = " << v_lhs_size << " and v_rhs = " << v_rhs_size
    << " in the operation v_lhs += op v_rhs do not match");
}

void DenseLinAlgPack::VopV_assert_sizes(size_type v_rhs1_size, size_type v_rhs2_size)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    v_rhs1_size != v_rhs2_size, std::length_error
    ,"VopV_assert_sizes(...) : The sizes of v_rhs1 and v_rhs2 "
    "in the operation v_rhs1 op v_rhs2 do not match");
}

void DenseLinAlgPack::Mp_M_assert_sizes(size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
  , size_type m_rhs_rows, size_type m_rhs_cols, BLAS_Cpp::Transp trans_rhs)
{
  if(		rows(m_lhs_rows,m_lhs_cols,trans_lhs) != rows(m_rhs_rows,m_rhs_cols,trans_rhs)
    ||	cols(m_lhs_rows,m_lhs_cols,trans_lhs) != cols(m_rhs_rows,m_rhs_cols,trans_rhs) )
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_M_assert_sizes(...) : The sizes of m_lhs and m_rhs "
      "in the operation op(m_lhs) += op op(m_rhs) do not match");
  }
}

void DenseLinAlgPack::MopM_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{
  if(		rows(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != rows(m_rhs2_rows,m_rhs2_cols,trans_rhs2)
    ||	cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != cols(m_rhs2_rows,m_rhs2_cols,trans_rhs2) )
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_M_assert_sizes(...) : The sizes of m_rhs1 and m_rhs2 "
      "in the operation op(m_rhs1) op op(m_rhs2) do not match");
  }
}

void DenseLinAlgPack::MtV_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols
  , BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{
  if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != v_rhs2_size)
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"MtV_assert_sizes(...) : The number of columns in "
      "m_rhs1 and the size of v_rhs2 in the operation v_lhs += op(m_rhs1) * v_rhs2 "
      "do not match");
}

void DenseLinAlgPack::Vp_MtV_assert_sizes(size_type v_lhs_size, size_type m_rhs1_rows
  , size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{
  if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != v_rhs2_size)
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Vp_MtV_assert_sizes(...) : The number of columns in"
      " m_rhs1 and the size of v_rhs2 in the operation v_lhs += op(m_rhs1) * v_rhs2"
      " do not match");
  if(rows(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != v_lhs_size)
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Vp_MtV_assert_sizes(...) : The number of rows in"
      " m_rhs1 and the size of v_lhs in the operation v_lhs += op(m_rhs1) * v_rhs2"
      " do not match");
}

void DenseLinAlgPack::MtM_assert_sizes(
    size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{
  if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != rows(m_rhs2_rows,m_rhs2_cols,trans_rhs2))
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"MtM_assert_sizes(...) : The number of columns in"
      " m_rhs1 and the number of rows in m_rhs2 in the operation"
      " op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
}

void DenseLinAlgPack::Mp_MtM_assert_sizes(
    size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
  , size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
  , size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{
  if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != rows(m_rhs2_rows,m_rhs2_cols,trans_rhs2))
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_MtM_assert_sizes(...) : The number of columns in"
      " m_rhs1 and the number of rows in m_rhs2 in the operation"
      " op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
  if(rows(m_lhs_rows,m_lhs_cols,trans_lhs) != rows(m_rhs1_rows,m_rhs1_cols,trans_rhs1))
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_MtM_assert_sizes(...) : The number of rows in"
      " m_lhs and the number of rows in m_rhs1 in the operation"
      " op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
  if(cols(m_lhs_rows,m_lhs_cols,trans_lhs) != cols(m_rhs2_rows,m_rhs2_cols,trans_rhs2))
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::length_error
      ,"Mp_MtM_assert_sizes(...) : The number of columns in"
      " m_lhs and the number of columns in m_rhs1 in the operation"
      " op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
}

#endif // LINALGPACK_CHECK_RHS_SIZES
