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

#include "RTOp_ROp_max_near_feas_step.h"
#include "RTOp_obj_value_vtbl.h"
#include "RTOp_obj_free_free.h"
#include "RTOp_get_reduct_op.hpp"

#include <stdlib.h>

/* */
/* Implementation functions */
/* */

/* Functions for the reduction target object */

static int get_targ_type_num_entries(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void* obj_data
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
  *num_values  = 2;
  *num_indexes = 0;
  *num_chars   = 0;
  return 0;
}

static int targ_obj_create(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* obj_data
  , RTOp_ReductTarget* targ_obj )
{
  const int mem_size = sizeof(struct RTOp_ROp_max_near_feas_step_reduct_obj_t);
  *targ_obj = malloc( mem_size );
  return vtbl->obj_reinit(vtbl,obj_data,*targ_obj);
}

static int targ_obj_reinit(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* obj_data
  , RTOp_ReductTarget targ_obj )
{
  struct RTOp_ROp_max_near_feas_step_reduct_obj_t
    *targ = (struct RTOp_ROp_max_near_feas_step_reduct_obj_t*)targ_obj;
  targ->alpha_pos = +1e+50; /* big enough? */
  targ->alpha_neg = -1e+50; /* big enough? */
  return 0;
}

static int targ_extract_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void *       obj_data
  ,void *             reduct_obj
  ,int                num_values
  ,RTOp_value_type    value_data[]
  ,int                num_indexes
  ,RTOp_index_type    index_data[]
  ,int                num_chars
  ,RTOp_char_type     char_data[]
  )
{
  struct RTOp_ROp_max_near_feas_step_reduct_obj_t *targ = NULL;
#ifdef RTOp_DEBUG
  assert( reduct_obj );
  assert( num_values  == 2 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
  assert( value_data );
#endif
  targ = (struct RTOp_ROp_max_near_feas_step_reduct_obj_t*)reduct_obj;
  value_data[0] = targ->alpha_pos;
  value_data[1] = targ->alpha_neg;
  return 0;
}

static int targ_load_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void*             obj_data
  ,int                     num_values
  ,const RTOp_value_type   value_data[]
  ,int                     num_indexes
  ,const RTOp_index_type   index_data[]
  ,int                     num_chars
  ,const RTOp_char_type    char_data[]
  ,void **                 reduct_obj
  )
{
  struct RTOp_ROp_max_near_feas_step_reduct_obj_t *targ = NULL;
#ifdef RTOp_DEBUG
  assert( *reduct_obj );
  assert( num_values  == 2 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
  assert( value_data );
#endif
  targ = (struct RTOp_ROp_max_near_feas_step_reduct_obj_t*)*reduct_obj;
  targ->alpha_pos = value_data[0];
  targ->alpha_neg = value_data[1];
  return 0;
}

static const struct RTOp_obj_type_vtbl_t targ_obj_vtbl =
{
  get_targ_type_num_entries
  ,targ_obj_create
  ,targ_obj_reinit
  ,RTOp_obj_free_free
  ,targ_extract_state
  ,targ_load_state
};

/* Other functions */

static int RTOp_ROp_max_near_feas_step_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  RTOp_value_type        beta;
  struct RTOp_ROp_max_near_feas_step_reduct_obj_t
    *targ = NULL;
  RTOp_index_type        sub_dim;
  const RTOp_value_type  *xl_val;
  ptrdiff_t              xl_val_s;
  const RTOp_value_type  *x_val;
  ptrdiff_t              x_val_s;
  const RTOp_value_type  *d_val;
  ptrdiff_t              d_val_s;
  const RTOp_value_type  *xu_val;
  ptrdiff_t              xu_val_s;
  register RTOp_index_type  k;
  RTOp_value_type  alpha;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 4 )
    return RTOp_ERR_INVALID_NUM_VECS;
  assert(vecs);
  if( num_targ_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  if( vecs[0].global_offset != vecs[1].global_offset
    || vecs[0].sub_dim != vecs[1].sub_dim
    || vecs[0].global_offset != vecs[2].global_offset
    || vecs[0].sub_dim != vecs[2].sub_dim
    || vecs[0].global_offset != vecs[3].global_offset
    || vecs[0].sub_dim != vecs[3].sub_dim
    )
    return RTOp_ERR_INCOMPATIBLE_VECS;

  /* */
  /* Get pointers to data */
  /* */

  /* beta */
  beta = *(RTOp_value_type*)obj_data;
  /* targ */
  targ = (struct RTOp_ROp_max_near_feas_step_reduct_obj_t*)targ_obj;
  /* sub_dim */
  sub_dim        = vecs[0].sub_dim;
  /* xl */
  xl_val         = vecs[0].values;
  xl_val_s       = vecs[0].values_stride;
  /* x */
  x_val          = vecs[1].values;
  x_val_s        = vecs[1].values_stride;
  /* d */
  d_val          = vecs[2].values;
  d_val_s        = vecs[2].values_stride;
  /* xu */
  xu_val         = vecs[3].values;
  xu_val_s       = vecs[3].values_stride;

  /* */
  /* If x has already been found to be infeasible just return */
  /* */
  if(targ->alpha_pos < 0.0)
    return 0; /* success! */

  /* */
  /* Perform the reduction operation. */
  /* */
  /* max alpha_pos, min alpha_neg s.t. xl - beta <= x + alpha*d <= xu + beta */
  /* */
  for( k = 0; k < sub_dim; ++k, xl_val += xl_val_s, x_val += x_val_s, d_val += d_val_s, xu_val += xu_val_s ) {
    if( *x_val < *xl_val - beta || *x_val > *xu_val + beta ) {
      targ->alpha_pos = -1.0; /* x is infeasible as is */
      return 0; /* success! */
    }
    if( *d_val != 0.0 ) {
      /* check lower bound */
      alpha = (*xl_val - beta - *x_val) / *d_val;
      if( ( alpha > 0.0 && alpha < targ->alpha_pos )
        || ( alpha == 0.0 && *d_val < 0.0 ) )
        targ->alpha_pos = alpha;
      if( ( alpha < 0.0 && -alpha < -targ->alpha_neg )
        || ( alpha == 0.0 && *d_val > 0.0 ) )
        targ->alpha_neg = alpha;
      /* check upper bound */
      alpha = (*xu_val + beta - *x_val) / *d_val;
      if( (alpha > 0.0 && alpha < targ->alpha_pos )
        || ( alpha == 0.0 && *d_val > 0.0 ) )
        targ->alpha_pos = alpha;
      if( ( alpha < 0.0 && -alpha < -targ->alpha_neg )
        || ( alpha == 0.0 && *d_val < 0.0 ) )
        targ->alpha_neg = alpha;
    }
  }

  return 0; /* success? */
}

static int reduce_reduct_objs(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data /* Can be NULL! */
  , RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj )
{
  const struct RTOp_ROp_max_near_feas_step_reduct_obj_t
    *i_targ = (const struct RTOp_ROp_max_near_feas_step_reduct_obj_t*)in_reduct_obj;
  struct RTOp_ROp_max_near_feas_step_reduct_obj_t
    *io_targ = (struct RTOp_ROp_max_near_feas_step_reduct_obj_t*)inout_reduct_obj;
  if(  i_targ->alpha_pos <  io_targ->alpha_pos )
    io_targ->alpha_pos = i_targ->alpha_pos;
  if( -i_targ->alpha_neg < -io_targ->alpha_neg )
    io_targ->alpha_neg = i_targ->alpha_neg;
  return 0;
}

INSERT_GET_REDUCT_OP_FUNCS(
  2,0,0,RTOp_ROp_max_near_feas_step_reduct_obj_t,reduce_reduct_objs
  ,targ_load_state,targ_extract_state
  ,external_reduct_op,get_reduct_op)

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_near_feas_step_vtbl =
{
  &RTOp_obj_value_vtbl /* use simple scalar value type for object instance data */
  ,&targ_obj_vtbl
  ,"ROp_max_near_feas_step"
  ,NULL
  ,RTOp_ROp_max_near_feas_step_apply_op
  ,reduce_reduct_objs
  ,get_reduct_op
};

/* Class specific functions */

int RTOp_ROp_max_near_feas_step_construct( RTOp_value_type beta, struct RTOp_RTOp* op )
{
  op->vtbl = &RTOp_ROp_max_near_feas_step_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  *((RTOp_value_type*)op->obj_data) = beta;
  return 0; /* success? */
}

int RTOp_ROp_max_near_feas_step_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->vtbl      = NULL;
  return 0; /* success? */
}

int RTOp_ROp_max_near_feas_step_set_beta( RTOp_value_type beta, struct RTOp_RTOp* op )
{
  *((RTOp_value_type*)op->obj_data) = beta;
  return 0; /* success? */
}

struct RTOp_ROp_max_near_feas_step_reduct_obj_t
RTOp_ROp_max_near_feas_step_val(RTOp_ReductTarget targ_obj)
{
  return *(struct RTOp_ROp_max_near_feas_step_reduct_obj_t*)targ_obj;
}
