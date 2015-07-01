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

#include "Teuchos_ConfigDefs.hpp"
#include "ExampleNLPDirectRTOps.h"
#include "RTOp_obj_null_vtbl.h"
#include "RTOp_obj_index_vtbl.h"

/* Implementation for RTOp_TOp_explnlp2_c_eval */

static int explnlp2_c_eval_apply_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, const int num_vecs, const struct RTOp_SubVector vecs[]
	, const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
	, RTOp_ReductTarget targ_obj )
{
	/* c */
	size_t                 sub_dim;
	RTOp_value_type        *c_val;
	ptrdiff_t              c_val_s;
	/* xD */
	const RTOp_value_type  *xD_val;
	ptrdiff_t              xD_val_s;
	/* xI */
	const RTOp_value_type  *xI_val;
	ptrdiff_t              xI_val_s;

	register RTOp_index_type  k;

	/*
   *Validate the input
   */
	if( num_vecs != 2 || vecs == NULL )
		return RTOp_ERR_INVALID_NUM_VECS;
	if( num_targ_vecs != 1 || targ_vecs == NULL )
		return RTOp_ERR_INVALID_NUM_TARG_VECS;
	if( targ_vecs[0].sub_dim != vecs[0].sub_dim
		|| targ_vecs[0].sub_dim != vecs[1].sub_dim
		|| targ_vecs[0].global_offset != vecs[0].global_offset
		|| targ_vecs[0].global_offset != vecs[1].global_offset )
		return RTOp_ERR_INCOMPATIBLE_VECS;

	/*
   * Get pointers to data
   */

	/* c */
	sub_dim       = targ_vecs[0].sub_dim;
	c_val         = targ_vecs[0].values;
	c_val_s       = targ_vecs[0].values_stride;
	/* xD */
	xD_val         = vecs[0].values;
	xD_val_s       = vecs[0].values_stride;
	/* xI */
	xI_val         = vecs[1].values;
	xI_val_s       = vecs[1].values_stride;

	/*
   * Compute c(j) = xI(i) * ( xD(i) - 1 ) - 10 * xD(i)
   */

	if( c_val_s == 1 && xD_val_s == 1 && xI_val_s == 1 ) {
		/* Slightly faster loop for unit stride vectors */
		for( k = 0; k < sub_dim; ++k, ++xI_val )
			*c_val++ = (*xD_val++) * (*xI_val - 1.0) - 10.0 * (*xI_val);
	}
	else {
		/* More general implementation for non-unit strides */
		for( k = 0; k < sub_dim; ++k, c_val+=c_val_s, xD_val+=xD_val_s, xI_val+=xI_val_s )
			*c_val = (*xD_val) * (*xI_val - 1.0) - 10.0 * (*xI_val);
	}

	return 0; /* success? */
}

const struct RTOp_RTOp_vtbl_t RTOp_TOp_explnlp2_c_eval_vtbl =
{
	&RTOp_obj_null_vtbl
	,&RTOp_obj_null_vtbl
	,"TOp_explnlp2_c_eval"
	,NULL
	,explnlp2_c_eval_apply_op
	,NULL
	,NULL
};

int RTOp_TOp_explnlp2_c_eval_construct( struct RTOp_RTOp* op )
{
	op->obj_data  = NULL;
	op->vtbl      = &RTOp_TOp_explnlp2_c_eval_vtbl;
	return 0;
}

int RTOp_TOp_explnlp2_c_eval_destroy( struct RTOp_RTOp* op )
{
	op->obj_data  = NULL;
	op->vtbl      = NULL;
	return 0;
}

/* Implementation for RTOp_TOp_explnlp2_calc_py_D */

static int explnlp2_calc_py_D_apply_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, const int num_vecs, const struct RTOp_SubVector vecs[]
	, const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
	, RTOp_ReductTarget targ_obj )
{
	size_t                 sub_dim;
	/* xD */
	const RTOp_value_type  *xD_val;
	ptrdiff_t              xD_val_s;
	/* xI */
	const RTOp_value_type  *xI_val;
	ptrdiff_t              xI_val_s;
	/* c */
	const RTOp_value_type  *c_val;
	ptrdiff_t              c_val_s;
	/* d */
	RTOp_value_type        *d_val;
	ptrdiff_t              d_val_s;
	/* py */
	RTOp_value_type        *py_val;
	ptrdiff_t              py_val_s;

	register RTOp_index_type  k;
	int                       all_unit_stride = 0;
	RTOp_value_type           denom;

	/* task */
	int task;
	assert(obj_data);
	task = *(int*)obj_data;
	assert(0 <= task && task <= 2);

	/*
   * Validate the input
   */
	if( ( (task == 0 || task == 1) && num_vecs != 2 )
		|| ( (task == 2) && num_vecs != 3 )
		|| vecs == NULL )
		return RTOp_ERR_INVALID_NUM_VECS;
	if( ( (task == 0 || task == 1) && num_targ_vecs != 1 )
		|| ( (task == 2) && num_targ_vecs != 2 )
		|| targ_vecs == NULL )
		return RTOp_ERR_INVALID_NUM_TARG_VECS;
	if( targ_vecs[0].sub_dim != vecs[0].sub_dim
		|| targ_vecs[0].sub_dim != vecs[1].sub_dim
		|| ( task == 2 && ( targ_vecs[0].sub_dim != vecs[2].sub_dim ) )
		|| ( task == 2 && ( targ_vecs[0].sub_dim != targ_vecs[1].sub_dim ) )
		|| targ_vecs[0].global_offset != vecs[0].global_offset
		|| targ_vecs[0].global_offset != vecs[1].global_offset
		|| ( task == 2 && (targ_vecs[0].global_offset != vecs[2].global_offset ) )
		|| ( task == 2 && ( targ_vecs[0].global_offset != targ_vecs[1].global_offset ) ) )
		return RTOp_ERR_INCOMPATIBLE_VECS;

	/*
   * Get pointers to data
   */

	sub_dim = vecs[0].sub_dim;

	k = 0;
	/* xD */
	xD_val         = vecs[k].values;
	xD_val_s       = vecs[k].values_stride;
	++k;
	if( task == 1 || task == 2 ) {
		/* xI */
		xI_val         = vecs[k].values;
		xI_val_s       = vecs[k].values_stride;
		++k;
	}
	if( task == 0 || task == 2 ) {
		/* c */
		c_val         = vecs[k].values;
		c_val_s       = vecs[k].values_stride;
		++k;
	}
	k = 0;
	if( task == 1 || task == 2 ) {
		/* d */
		d_val         = targ_vecs[k].values;
		d_val_s       = targ_vecs[k].values_stride;
		++k;
	}
	if( task == 0 || task == 2 ) {
		/* py */
		py_val         = targ_vecs[k].values;
		py_val_s       = targ_vecs[k].values_stride;
		++k;
	}

	/* Determine if all the vectors have unit stride! */
	all_unit_stride = 1;
	for( k = 0; k < num_vecs && !all_unit_stride; ++k )
		if( vecs[k].values_stride != 1 )
			all_unit_stride = 0;
	for( k = 0; k < num_targ_vecs && !all_unit_stride; ++k )
		if( targ_vecs[k].values_stride != 1 )
			all_unit_stride = 0;

	/*
   * Compute py and/or D
   */

	if( all_unit_stride) {
		if(task == 0) {
			/* Compute py only */
			for( k = 0; k < sub_dim; ++k )
				*py_val++ = *c_val++ / ( 1.0 - *xI_val++ );
		}
		if(task == 1) {
			/* Compute D only */
			for( k = 0; k < sub_dim; ++k )
				*d_val++ = ( *xD_val++ - 10.0 ) / ( 1.0 - *xI_val++ );
		}
		if(task == 2) {
			/* Compute py and D */
			for( k = 0; k < sub_dim; ++k ) {
				denom = ( 1.0 - *xI_val++ );
				*d_val++ = ( *xD_val++ - 10.0 ) / denom;
				*py_val++ = *c_val++ / denom;
			}
		}
	}
	else {
		assert(0); /* ToDo: Implement if needed! */
	}

	return 0;
}

const struct RTOp_RTOp_vtbl_t RTOp_TOp_explnlp2_calc_py_D_vtbl =
{
	&RTOp_obj_index_vtbl
	,&RTOp_obj_null_vtbl
	,"TOp_explnlp2_calc_py_D"
	,NULL
	,explnlp2_calc_py_D_apply_op
	,NULL
	,NULL
};

int RTOp_TOp_explnlp2_calc_py_D_construct( int task, struct RTOp_RTOp* op )
{
	int result;
#ifdef RTOp_DEBUG
	assert( 0 <= task && task <= 2 );
#endif	
	op->obj_data  = NULL;
	op->vtbl      = &RTOp_TOp_explnlp2_calc_py_D_vtbl;
	result = op->vtbl->obj_data_vtbl->obj_create( NULL, NULL, &op->obj_data );
	if(result != 0) return result;
#ifdef RTOp_DEBUG
	assert(op->obj_data);
#endif
	*((int*)op->obj_data) = task;
	return 0;
}

int RTOp_TOp_explnlp2_calc_py_D_set_task( int task, struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
	assert( 0 <= task && task <= 2 );
	assert(op->obj_data);
#endif
	*((int*)op->obj_data) = task;
	return 0;
}

int RTOp_TOp_explnlp2_calc_py_D_destroy( struct RTOp_RTOp* op )
{
	int result;
	result = op->vtbl->reduct_vtbl->obj_free( NULL, NULL, &op->obj_data );
	if(result != 0) return result;
	op->obj_data  = NULL;
	op->vtbl      = NULL;
	return 0;
}
