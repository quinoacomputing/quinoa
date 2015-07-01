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

#include "RTOp_SparseSubVector.h"

void RTOp_sparse_sub_vector(
	RTOp_index_type global_offset, RTOp_index_type sub_dim
	,RTOp_index_type sub_nz, const RTOp_value_type values[], ptrdiff_t values_stride
	,const RTOp_index_type indices[], ptrdiff_t indices_stride
	,ptrdiff_t local_offset, int is_sorted
	,struct RTOp_SparseSubVector *sub_vec
	)
{
	/* Validate input */
#ifdef RTOp_DEBUG
	assert( sub_vec );
	assert(
		( sub_nz && ( values != NULL && indices != NULL && indices_stride != 0 && sub_nz <= sub_dim ) )
		|| !sub_nz || ( sub_nz == sub_dim && indices == NULL )
		);
#endif
	/* Set members */
	sub_vec->global_offset  = global_offset;
	sub_vec->sub_dim        = sub_dim;
	sub_vec->sub_nz         = sub_nz;
	sub_vec->values         = values;
	sub_vec->values_stride  = values_stride;
	sub_vec->indices        = indices;
	sub_vec->indices_stride = indices_stride;
	sub_vec->local_offset   = local_offset;
	sub_vec->is_sorted      = is_sorted;
}

void RTOp_sparse_sub_vector_null( struct RTOp_SparseSubVector *sub_vec )
{
	sub_vec->global_offset  = 0;
	sub_vec->sub_dim        = 0;
	sub_vec->sub_nz         = 0;
	sub_vec->values         = NULL;
	sub_vec->values_stride  = 0;
	sub_vec->indices        = NULL;
	sub_vec->indices_stride = 0;
	sub_vec->local_offset   = 0;
	sub_vec->is_sorted      = 0;
}

void RTOp_sparse_sub_vector_from_dense(
	const struct RTOp_SubVector     *sub_vec
	,struct RTOp_SparseSubVector    *spc_sub_vec
	)
{
	spc_sub_vec->global_offset  = sub_vec->global_offset;
	spc_sub_vec->sub_dim        = sub_vec->sub_dim;
	spc_sub_vec->sub_nz         = sub_vec->sub_dim;
	spc_sub_vec->values         = sub_vec->values;
	spc_sub_vec->values_stride  = sub_vec->values_stride;
	spc_sub_vec->indices        = NULL;
	spc_sub_vec->indices_stride = 0;
	spc_sub_vec->local_offset   = 0;
	spc_sub_vec->is_sorted      = 0;
}
