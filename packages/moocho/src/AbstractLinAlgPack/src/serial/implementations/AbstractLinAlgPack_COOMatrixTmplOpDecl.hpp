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

#ifndef COO_MATRIX_TMPL_OP_DECL_H
#define COO_MATRIX_TMPL_OP_DECL_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \defgroup COOMatrixTmplOp_grp Linear algebra operations for coordinate matrices.
 *
 * This is a basic set of BLAS like linear algebra operations with sparse
 * coordinate (COO) matrices accesses through a template interface
 * called COOMatrixTemplateInterface.  Through this template interface the client
 * can inquire as to the dimensions of the matrix (#rows()#, #cols()#) and how many
 * nonzero elements it has (#nz()#).  The nonzero elements are accessed using forward
 * iterators returned from #begin()# and #end()#.  The iterator must return a type
 * with an interface compatable with the SparseCOOElementTemplateInterface specification.
 *
 * The specifications for these template interfaces are given below:\\
 *
 \begin{verbatim}
  class SparseCOOElementTemplateInterface {
  public:
    typedef ...		value_type;
    typedef ...		index_type;

    value_type& value();
    value_type value() const;
    index_type row_i() const;
    index_type col_j() const;
  };

  class COOMatrixTemplateInterface {
  public:
    typedef ...		size_type;
    typedef ...		difference_type;
    typedef ...		element_type;		// SparseCOOElementTemplateInterface compliant
    typedef ...		iterator;			// returns an element_type
    typedef ...		const_iterator;		// returns a const element_type

    size_type rows() const;
    size_type cols() const;
    size_type nz() const;
    difference_type row_offset() const;
    difference_type col_offset() const;
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;
  };
 \end{verbatim}
 *
 * The nonzero elements are specified by the triplet (val,i,j) where:\\
 *
 * #val[k]	= (begin() + k)->value()#\\
 * #i[k]	= (begin() + k)->row_i() + row_offset()#\\
 * #j[k]	= (begin() + k)->col_j() + col_offset()#\\
 *
 * for #k# = #0#, ..., #nz() - 1#.
 *
 * The iterator returned by the begin() functions needs to be a forward
 * iterator that supports <tt>++itr</tt>.  These elements need not be sorted but
 * the operations may perform better if they are.  What order will dependent on
 * which operation and which options.  See the individual operations to determine
 * which sorting orders are best.
 *
 * Actually these linear algebra functions only use the constant interface to
 * these tempalte interfaces so the iterators only need to support the constant
 * interface.
 *
 * These functions use the same basic naming sceme as DenseLinAlgPack with a
 * few exceptions.  All of the operations are += instead of = because of the
 * nature of the sparse matrix calculates so instead of #V_# and #M_# starting
 * out the function names they begin with #Vp_# and #Mp_# with #p_# used
 * to signify += instead of just =.  Also #COOM# replaces #M# where a 
 * COO matrrix replaces the dense DMatrix.  This is because the overload
 * resolution rules would mess up the implicit type conversions that
 * go on in DenseLinAlgPack if the exact same names where used.
 */
//@{

/// gms_lhs += alpha * op(coom_rhs) (time = O(coom_rhs.nz(), space = O(1))
template<class T_COOM>
void Mp_StCOOM(DMatrixSlice* gms_lhs, value_type alpha, const T_COOM& coom_rhs
  , BLAS_Cpp::Transp trans_rhs);

/// vs_lhs += alpha * op(coom_rhs1) * vs_rhs2 (BLAS xGEMV) (time = O(coom_rhs.nz(), space = O(1))
template<class T_COOM>
void Vp_StCOOMtV(DVectorSlice* vs_lhs, value_type alpha, const T_COOM& coom_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2);

/// gms_lhs += alpha * op(coom_rhs1) * op(gms_rhs2) (right) (BLAS xGEMM)
template<class T_COOM>
void Mp_StCOOMtM(DMatrixSlice* gms_lhs, value_type alpha, const T_COOM& coom_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2);

/// gms_lhs += alpha * op(gms_rhs1) * op(coom_rhs2) (left) (BLAS xGEMM)
template<class T_COOM>
void Mp_StMtCOOM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSlice& gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const T_COOM& coom_rhs2, BLAS_Cpp::Transp trans_rhs2);

// / gms_lhs = alpha * op(coom_rhs1) * op(sym_rhs2) (right) (BLAS xSYMM)
//template<class T_COOM>
//void Mp_StCOOMtSM(DMatrixSlice* gms_lhs, value_type alpha, const T_COOM& coom_rhs1
//	, BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceSym& sym_rhs2, BLAS_Cpp::Transp trans_rhs2);

// / gms_lhs = alpha * op(sym_rhs1) * op(coom_rhs2) (left) (BLAS xSYMM)
//template<class T_COOM>
//void Mp_StSMtCOOM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceSym& sym_rhs1
//	, BLAS_Cpp::Transp trans_rhs1, const T_COOM& coom_rhs2, BLAS_Cpp::Transp trans_rhs2);

// / gms_lhs = alpha * op(coom_rhs1) * op(tri_rhs2) (right) (BLAS xTRMM)
//template<class T_COOM>
//void Mp_StCOOMtSM(DMatrixSlice* gms_lhs, value_type alpha, const T_COOM& coom_rhs1
//	, BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceTri& tri_rhs2, BLAS_Cpp::Transp trans_rhs2);

// / gms_lhs = alpha * op(tri_rhs1) * op(coom_rhs2) (left) (BLAS xTRMM)
//template<class T_COOM>
//void Mp_StSMtCOOM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceTri& tri_rhs1
//	, BLAS_Cpp::Transp trans_rhs1, const T_COOM& coom_rhs2, BLAS_Cpp::Transp trans_rhs2);

//@}

} // end namespace AbstractLinAlgPack

#endif	// COO_MATRIX_TMPL_OP_DECL_H
