/*
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
*/

#ifndef RTOP_SPARSE_SUB_VECTOR_H
#define RTOP_SPARSE_SUB_VECTOR_H

#include <stddef.h>

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* */
/** Struct for a (sparse or dense) sub-vector.
 *
 *	Sparse and dense local vectors are supported as follows:
  *
  *	A dense vector <tt>vec</tt> is identified by <tt>vec.sub_dim == vec.sub_nz</tt>
  * and <tt>vec.indices == NULL</tt> in which case
  *	<tt>vec.indices_stride</tt>, <tt>vec.local_offset</tt> and <tt>vec.is_sorted</tt>
  * are ignored.  For a dense sub-vector <tt>vec</tt>, the corresponding entries
 *	in the global vector <tt>x(j)</tt> (one based) are as follows:
 \verbatim

	x( vec.global_offset + k )
		= vec.values[ vec.value_stride * (k-1) ]

	for k = 1,...,vec.sub_dim
 \endverbatim
 * The stride member <tt>vec.value_stride</tt> may be positive (>0), negative (<0)
 * or even zero (0).  A negative stride <tt>vec.value_stride < 0</tt> allows a
 * reverse traversal of the elements in <tt>vec.values[]</tt>.  A zero stride
 * <tt>vec.value_stride == 0</tt> allows a vector with all the elements the same.
 *
 *	A sparse vector is identified by <tt>vec.sub_dim > vec.sub_nz</tt>
 * or <tt>vec.indices != NULL</tt>
 * in which case all the fields in the structure are meaningful.
 *	The corresponding elements in the global vector <tt>x(j)</tt>
 * defined as:
 \verbatim

	x( vec.global_offset + vec.local_offset + vec.indices[vec.indices_stride*(k-1)] )
		= vec.values[vec.value_stride*(k-1)]

	for k = 1,...,vec.sub_nz
 \endverbatim
 * If <tt>vec.sub_nz == 0</tt> then it is allowed for <tt>vec.indices == NULL</tt>.
 * If <tt>vec.sub_dim > vec.sub_nz > 0</tt> then <tt>vec.indices != NULL</tt> must be true.
 *
  * A sparse sub-vector may be sorted (<tt>vec.is_sorted!=0</tt>) or
  * unsorted (<tt>vec.is_sorted==0</tt>) but the indices <tt>vec.indices[k]</tt>
  * must be unique.  A sorted vector (<tt>vec.is_sorted!=0</tt>) means that
  * the indices are in ascending order:
  \verbatim

	vec.indices[vec.indices_stride*(k-1)] < vec.indices[vec.indices_stride*(k)]

	for k = 1,...,vec.sub_nz-1
 \endverbatim
 * The member <tt>vec.local_offset</tt> is used to shift the values in <tt>vec.indices[]</tt>
 * to be in range of the local sub-vector.  In other words:
 \verbatim
	
	1 <= vec.local_offset + vec.indices[vec.indices_stride*(k-1)] <= vec.sub_nz

	for k = 1...vec.sub_nz
 \endverbatim
 * The member <tt>vec.value_stride</tt> may be positive (>0), negative (<0) or zero (0).
 * However, the member <tt>vec.indices_stride</tt> may be may be positive (>0)
 * or negative (<0) but not zero (0).  Allowing <tt>vec.indices_stride == 0</tt>
 * would mean that a vector would have <tt>vec.sub_nz</tt> nonzero elements with
 * all the same value and all the same indexes and non-unique indices are
 * not allowed.  Allowing non-unique indexes would make some operations
 * (e.g. dot product) very difficult to implement and therefore can not
 * be allowed.  A sparse vector where <tt>vec.value_stride == 0</tt> is one
 * where all of the nonzeros have the value <tt>vec.values[0]</tt>.  If
 * <tt>vec.sub_nz == 0</tt> for a sparse vector then it is allowed for
 * <tt>vec.values == NULL</tt> and <tt>vec.indices == NULL</tt>.
 *
 *	This specification allows a lot of flexibility in determining
 * how the vectors are laid out in memory.  However, allowing vectors to be
 * sparse and unsorted may make many user defined operations
 * considerably harder and expensive to implement.
 *
 * To avoid making mistakes in setting the members of this struct use
 * one of the helper functions <tt>RTOp_sparse_sub_vector_from_dense()</tt>,
 * <tt>RTOp_sparse_sub_vector()</tt> or <tt>RTOp_sub_vector_null()</tt>.
 */
struct RTOp_SparseSubVector {
	/* Offset for the sub-vector into the global vector */
	RTOp_index_type                  global_offset;
	/* Dimension of the sub-vector */
	RTOp_index_type                  sub_dim;
	/* Number of nonzero elements (<tt>sub_nz == sub_dim</tt> for dense vectors) */
	RTOp_index_type                  sub_nz;
	/* Array (size min{|<tt>value_stride*sub_nz</tt>|,1}) for the values in the vector */
	const RTOp_value_type            *values;
	/* Stride between elements in <tt>values[]</tt> */
	ptrdiff_t                        values_stride;
	/* */
	/** Array (size min{|<tt>indices_stride*sub_nz</tt>|,1} if not <tt>NULL</tt>) for the
	  * indices of the nonzero elements in the vector (sparse vectors only)
	  */
	const RTOp_index_type            *indices;
	/* Stride between indices in indices[] (sparse vectors only) */
	ptrdiff_t                        indices_stride;
	/* Offset of indices[] into local sub-vector (sparse vectors only) */
	ptrdiff_t                        local_offset;
	/* If <tt>is_sorted == 0</tt> then the vector is not sorted, otherwise it is sorted (sparse vectors only) */
	int                              is_sorted;
};

/* */
/** Set the members for sparse sub-vector.
  */
void RTOp_sparse_sub_vector(
	RTOp_index_type global_offset, RTOp_index_type sub_dim, RTOp_index_type sub_nz
	,const RTOp_value_type values[], ptrdiff_t values_stride
	,const RTOp_index_type indices[], ptrdiff_t indices_stride
	,ptrdiff_t local_offset, int is_sorted
	,struct RTOp_SparseSubVector *sub_vec
	);
/* */
/** Initialize a sub-vector argument to null.
  */
void RTOp_sparse_sub_vector_null( struct RTOp_SparseSubVector *sub_vec );

/* */
/** Copy the elements from a RTOp_SubVector to a RTOp_SparseSubVector.
  */
void RTOp_sparse_sub_vector_from_dense(
	const struct RTOp_SubVector     *sub_vec
	,struct RTOp_SparseSubVector    *spc_sub_vec
	);

#ifdef __cplusplus
}
#endif

#endif /* RTOP_SPARSE_SUB_VECTOR_H */
