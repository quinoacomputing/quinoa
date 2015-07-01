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

#include <math.h>

#include "RTOp_ROp_max_abs_ele.h"
#include "RTOp_obj_null_vtbl.h"
#include "RTOp_obj_value_vtbl.h"
#include "RTOp_obj_free_free.h"
#include "RTOp_get_reduct_op.hpp"

/* */
/* Implementation functions */
/* */

/* Selected functions that are used to implement exteral_reduct_op */

static int CALL_API targ_extract_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void *       instance_data
  ,void *             obj
  ,int                num_values
  ,RTOp_value_type    value_data[]
  ,int                num_indexes
  ,RTOp_index_type    index_data[]
  ,int                num_chars
  ,RTOp_char_type     char_data[]
  )
{
  struct RTOp_value_index_type* vi_obj;
#ifdef RTOp_DEBUG
  assert( obj );
  assert( num_values  == 1 );
  assert( num_indexes == 1 );
  assert( num_chars   == 0 );
#endif
  vi_obj = (struct RTOp_value_index_type*)obj;
  value_data[0] = vi_obj->value;
  index_data[0] = vi_obj->index;
  return 0;
}

static int CALL_API targ_load_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void *            instance_data
  ,int                     num_values
  ,const RTOp_value_type   value_data[]
  ,int                     num_indexes
  ,const RTOp_index_type   index_data[]
  ,int                     num_chars
  ,const RTOp_char_type    char_data[]
  ,void **                 obj
  )
{
  struct RTOp_value_index_type* vi_obj;
#ifdef RTOp_DEBUG
  assert( obj );
  assert( *obj );
  assert( num_values  == 1 );
  assert( num_indexes == 1 );
  assert( num_chars   == 0 );
#endif
  vi_obj = (struct RTOp_value_index_type*)*obj;
  vi_obj->value = value_data[0];
  vi_obj->index = index_data[0];
  return 0;
}

/* Other functions */

static int RTOp_ROp_max_abs_ele_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  ,const int num_vecs, const struct RTOp_SubVector vecs[]
  ,const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  ,RTOp_ReductTarget targ_obj
  )
{
  /* */
  /* Declare local variables */
  /* */

  /* targ */
  struct RTOp_value_index_type
    *targ = NULL;
  /* global_off */
  size_t                 global_offset;
  /* sub_dim */
  size_t                 sub_dim;
  /* v */
  const RTOp_value_type  *v_val = NULL;
  ptrdiff_t              v_val_s;

  register size_t  k;
  RTOp_index_type  i;
  RTOp_value_type  abs_v_i;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 1 )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;

  /* */
  /* Get pointers to the data */
  /* */

  /* targ */
  targ            = (struct RTOp_value_index_type*)targ_obj;
  /* global_off */
  global_offset   = vecs[0].global_offset;
  /* sub_dim */
  sub_dim         = vecs[0].sub_dim;
  /* v */
  v_val           = vecs[0].values;
  v_val_s         = vecs[0].values_stride;

  /* */
  /* Perform the reduction operation. */
  /* */

  i = global_offset + 1;
  for( k = 0; k < sub_dim; ++k, ++i, v_val += v_val_s ) {
    abs_v_i = fabs(*v_val);
    if( abs_v_i > targ->value || ( abs_v_i == targ->value && i < targ->index ) || targ->index == 0 ) {
      targ->value = *v_val;
      targ->index = i;
    }
  }

  return 0; /* success? */
}

static int reduce_reduct_objs(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data /* Can be NULL! */
  , RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj )
{
  const struct RTOp_value_index_type
    *i_targ = (const struct RTOp_value_index_type*)in_reduct_obj;
  struct RTOp_value_index_type
    *io_targ = (struct RTOp_value_index_type*)inout_reduct_obj;
  RTOp_value_type
    i_abs_val  = fabs(i_targ->value),
    io_abs_val = fabs(io_targ->value);
  if(
    ( i_abs_val > io_abs_val )
    ||
    ( i_abs_val > io_abs_val && i_targ->index < io_targ->index )
    )
  {
    io_targ->value     = i_targ->value;
    io_targ->index     = i_targ->index;
  }
  return 0;
}

INSERT_GET_REDUCT_OP_FUNCS(
  1,1,0,RTOp_value_index_type,reduce_reduct_objs
  ,targ_load_state,targ_extract_state
  ,external_reduct_op,get_reduct_op)

const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_abs_ele_vtbl =
{
  &RTOp_obj_null_vtbl
  ,&RTOp_obj_value_index_vtbl
  ,"ROp_max_abs_ele"
  ,NULL
  ,RTOp_ROp_max_abs_ele_apply_op
  ,reduce_reduct_objs
  ,get_reduct_op
};

/* Class specific functions */

int RTOp_ROp_max_abs_ele_construct( struct RTOp_RTOp* op )
{
  op->vtbl     = &RTOp_ROp_max_abs_ele_vtbl;
  op->obj_data = NULL;
  return 0; /* success? */
}

int RTOp_ROp_max_abs_ele_destroy( struct RTOp_RTOp* op )
{
  op->vtbl     = NULL;
  op->obj_data = NULL;
  return 0; /* success? */
}

struct RTOp_value_index_type
RTOp_ROp_max_abs_ele_val(RTOp_ReductTarget targ_obj)
{
  return *(struct RTOp_value_index_type*)targ_obj;
}
