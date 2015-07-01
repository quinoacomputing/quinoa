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

#ifndef MATRIX_VECTOR_TEMPLATE_OP_H
#define MATRIX_VECTOR_TEMPLATE_OP_H

#include <stdexcept>
#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** @name {\bf Templated Matrix-DVector Operations}.
  *
  * These are the declarations (AbstractLinAlgPack_MatrixVectorTemplateOp.hpp) for template functions for
  * performing selected matrix-vector and matrix-matrix operations.  The templated
  * matrix type must have the following interface:
  * \begin{itemize}
  * \item #T_Matrix::size_type# - Member type for the type returned by #rows# and #cols#
  * \item #T_Matrix::row_type# - Member type for the type returned by #row(i)#
  * \item #T_Matrix::col_type# - Member type for the type returned by #col(j)#
  * \item #T_Matrix::size_type T_Matrix::rows() const#
  * \item #T_Matrix::size_type T_Matrix::cols() const#
  * \item #const T_Matrix::row_type& T_Matrix::row(size_type i) const# \
  *		- Returns a reference to the ith row (1,2,...,#T_Matrix::rows()#).
  * \item #const T_Matrix::col_type& T_Matrix::col(size_type j) const# \ 
  *		- Returns a reference to the ith row (1,2,...,#T_Rect_Matrix::rows()#).
  * \end{itemize}
  *
  * In addition, the following operations (UML notation) must be accessible for the types
  * #T_Matrix::row_type# and #T_Matrix::col_type# when these template functions
  * are instantiated:
  *
  * assign(v_lhs:DVector&,v_rhs:T_Matrix::const row_type&)
  * assign(vs_lhs:DVectorSlice&,v_rhs:const T_Matrix::row_type&)
  * dot(vs_rhs1:const DVectorSlice&, vs_rhs2: const T_Matrix::row_type&):value_type 
  * dot(vs_rhs1:const DVectorSlice&, vs_rhs2: const T_Matrix::col_type&):value_type 
  */
//@{

// ////////////////////////
// assign to a dense matrix

/// gm_lhs = T_M (templated matrix type T_M)
template<class T_Matrix>
void assign(DMatrix& gm_lhs, const T_Matrix& gm_rhs, BLAS_Cpp::Transp trans_rhs);

/// gms_lhs = T_M (templated matrix type T_M)
template<class T_Matrix>
void assign(DMatrixSlice& gms_lhs, const T_Matrix& gm_rhs, BLAS_Cpp::Transp trans_rhs);

// /////////////////////////////
// Matrix-DVector multiplication

/// v_lhs = T_M * vs_lhs (templated matrix type T_M)
template<class T_Matrix>
void V_MtV(DVector& v_lhs, const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& vs_rhs2);

/// vs_lhs = T_M * vs_lhs (templated matrix type T_M)
template<class T_Matrix>
void V_MtV(DVectorSlice& v_lhs, const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& vs_rhs2);

//@}

}

#endif // MATRIX_VECTOR_TEMPLATE_OP_H
