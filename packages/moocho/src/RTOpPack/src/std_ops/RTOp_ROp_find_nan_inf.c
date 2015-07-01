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

#include "RTOp_ROp_find_nan_inf.h"
#include "RTOp_obj_null_vtbl.h"
#include "RTOp_obj_value_index_vtbl.h"
#include "check_nan_inf.h"
#include "RTOp_get_reduct_op.hpp"

/* Implementation functions */

static int RTOp_ROp_find_nan_inf_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  struct RTOp_value_index_type
                         *targ;
  RTOp_index_type        global_offset;
  RTOp_index_type        sub_dim;
  const RTOp_value_type  *v0_val;
  ptrdiff_t              v0_val_s;
  RTOp_index_type        i;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 1 )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  assert(targ_obj);
  assert(vecs);

  /* */
  /* Get pointers to data */
  /* */

  /* {value,index} */
  targ = (struct RTOp_value_index_type*)targ_obj;

  /* v0 */
  global_offset  = vecs[0].global_offset;
  sub_dim        = vecs[0].sub_dim;
  v0_val         = vecs[0].values;
  v0_val_s       = vecs[0].values_stride;

  /* */
  /* Find a NaN or Inf element! */
  /* */

  for( i = global_offset + 1; i <= global_offset + sub_dim; ++i, v0_val += v0_val_s ) {
    if( RTOp_is_nan_inf(*v0_val) && ( targ->index == 0 || i < targ->index ) ) {
      targ->value = *v0_val;
      targ->index = i;
      break;
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
  if(  io_targ->index == 0 || ( i_targ->index != 0 && i_targ->index < io_targ->index ) ) {
    io_targ->index = i_targ->index;
    io_targ->value = i_targ->value;
  }
  return 0;
}

INSERT_GET_REDUCT_OP_FUNCS(
  1,1,0,RTOp_value_index_type,reduce_reduct_objs
  ,RTOp_obj_value_index_vtbl.load_state
  ,RTOp_obj_value_index_vtbl.extract_state
  ,external_reduct_op,get_reduct_op)

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_ROp_find_nan_inf_vtbl =
{
  &RTOp_obj_null_vtbl         /* use null type for instance data */
  ,&RTOp_obj_value_index_vtbl /* use {value,index} type for target object */
  ,"ROp_find_nan_inf"
  ,NULL
  ,RTOp_ROp_find_nan_inf_apply_op
  ,reduce_reduct_objs
  ,get_reduct_op
};

/* Class specific functions */

int RTOp_ROp_find_nan_inf_construct( struct RTOp_RTOp* op )
{
  op->obj_data = NULL;
  op->vtbl     = &RTOp_ROp_find_nan_inf_vtbl;
  return 0;
}

int RTOp_ROp_find_nan_inf_destroy( struct RTOp_RTOp* op )
{
  op->obj_data = NULL;
  op->vtbl     = NULL;
  return 0;
}

struct RTOp_ROp_find_nan_inf_reduct_obj_t
RTOp_ROp_find_nan_inf_val(RTOp_ReductTarget targ_obj)
{
  const struct RTOp_value_index_type                *value_index;
  struct RTOp_ROp_find_nan_inf_reduct_obj_t         return_obj;
#ifdef RTOp_DEBUG
  assert(targ_obj != RTOp_REDUCT_OBJ_NULL);
#endif
  value_index = (const struct RTOp_value_index_type*)targ_obj;
  return_obj.v0_i = value_index->value;
  return_obj.i    = value_index->index;
  return return_obj;
}
