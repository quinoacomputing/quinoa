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
// Declarations for sparse vector operations.
// 

#ifndef SPARSE_VECTOR_OP_DECL_H
#define SPARSE_VECTOR_OP_DECL_H

#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"

namespace AbstractLinAlgPack {

/** @name Template operations for sparse vectors.
  *
  * These functions implement level 1 and 2 BLAS like
  * linear algebra operations on unsorted sparse
  * vectors.  These templated sparse vector objects
  * give information about the sparse vector (size, nonzeros etc.)
  * and give access to the sparse elements as iterators.
  * The iterators yield sparse elements that give the elements
  * indice in the full vector and its value.
  *
  * The specification for these interfaces is as follows:
  *
  \begin{verbatim}
  class SparseElementTemplateInterface {
  public:
    typedef ....	value_type;
    typedef ....	indice_type;

    value_type&		value();
    value_type		value() const;
    indice_type		indice() const;
  };

  class SparseVectorTemplateInterface {
  public:
    typedef ...		difference_type;
    typedef ...		element_type;	// SparseElementTemplateInterface compliant
    typedef ...		iterator;		// *(iter) yields a element_type
    typedef ...		const_iterator;	// *(iter) yields a const element_type
    typedef ...		reverse_iterator;	// *(iter) yields a element_type
    typedef ...		const_reverse_iterator;	// *(iter) yields a const element_type

    // Information
    size_type		size() const;	// size of the full vector
    size_type		nz() const;		// number of nonzero elements
    difference_type	offset() const;	// ith real real indice = begin()[i-1] + offset()
    bool			is_sorted() const;	// true if elements are sorted by indice

    // iterate forward (sorted) through elemlents
    iterator		begin();
    const_iterator	begin() const;
    iterator		end();
    const_iterator	end() const;

    // iterate backward (sorted) through elemlents
    reverse_iterator		rbegin();
    const_reverse_iterator	rbegin() const;
    reverse_iterator		rend();
    const_reverse_iterator	rend() const;
  };
/end{verbatim}
  *
  * In all of these functions where we have some operation that yields a dense vector
  * being added to another function such as:
  *
  * v_lhs = operation + vs_rhs2
  *
  * it is allowed that v_lhs.overlap(vs_rhs2) == SAME_MEM.  In this case no unnecesary
  * operations will be performed.  Also, it is up the the user to ensure that
  * there is not an alias problem where if v_lhs is the same as vs_rhs2 and vs_rhs2
  * is also used in the operation.  This has undefined results.  If a future version
  * of the library this may be handeled but for now it is not.
  *
  * These operations use the same nameing convensions as those for DVector and
  * DVectorSlice in DenseLinAlgPack with the acception that the sparse vectors
  * are given the tag #SV# instead of #V# so as to not destroy the intended
  * behavior of the operations in DenseLinAlgPack and the implicit conversion of
  * a DVector to a DVectorSlice.
  *
  * It should be noted that these operations will be more efficient for large
  * dense vectors and sparse vectors with many nonzero elements if the sparse
  * elements are sorted by indice.  This is because in many operations the elements of the
  * dense vectors are accessed using random access and this could cause
  * virtual page thrashing if the nonzero sparse elements are not sorted.
  */
//@{

/// result = dot(vs_rhs1,sv_rhs2) (BLAS xDOT)
template<class T_SpVec>
value_type dot_V_SV(const DVectorSlice& vs_rhs1, const T_SpVec& sv_rhs2);

/// result = dot(sv_rhs1,vs_rhs2) (BLAS xDOT)
template<class T_SpVec>
value_type dot_SV_V(const T_SpVec& sv_rhs1, const DVectorSlice& vs_rhs2);

/// result = ||sv_rhs||1 (BLAS xASUM)
template<class T_SpVec>
value_type norm_1_SV(const T_SpVec& sv_rhs);

/// result = ||sv_rhs||2 (BLAS xNRM2)
template<class T_SpVec>
value_type norm_2_SV(const T_SpVec& sv_rhs);

/// result = ||sv_rhs||inf (BLAS IxAMAX)
template<class T_SpVec>
value_type norm_inf_SV(const T_SpVec& sv_rhs);

/// result = max(sv_rhs)
template<class T_SpVec>
value_type max_SV(const T_SpVec& sv_rhs);

/// result = min(sv_rhs)
template<class T_SpVec>
value_type min_SV(const T_SpVec& sv_rhs);

/// sv_lhs *= alpha (BLAS xSCAL)
template<class T_SpVec>
void Vt_S( T_SpVec* sv_lhs, value_type alpha );

/// vs_lhs += alpha * sv_rhs (BLAS xAXPY)
template<class T_SpVec>
void Vp_StSV(DVectorSlice* vs_lhs, value_type alpha, const T_SpVec& sv_rhs);

/// vs_lhs += alpha * op(gms_rhs1) * sv_rhs2 (BLAS xGEMV)
template<class T_SpVec>
void Vp_StMtSV(DVectorSlice* vs_lhs, value_type alpha, const DMatrixSlice& gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const T_SpVec& sv_rhs2);

/// vs_lhs += alpha * op(tri_gms_rhs1) * sv_rhs2 (BLAS xTRMV)
template<class T_SpVec>
void Vp_StMtSV(DVectorSlice* vs_lhs, value_type alpha, const DMatrixSliceTri& tri_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const T_SpVec& sv_rhs2);

/// vs_lhs += alpha * op(sym_gms_rhs1) * sv_rhs2 (BLAS xSYMV)
template<class T_SpVec>
void Vp_StMtSV(DVectorSlice* vs_lhs, value_type alpha, const DMatrixSliceSym& sym_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const T_SpVec& sv_rhs2);

//@}

} // end namespace AbstractLinAlgPack

#endif // SPARSE_VECTOR_OP_DECL_H
