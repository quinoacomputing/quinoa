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

#include "RTOp_TOp_force_in_bounds.h"
#include "RTOp_obj_null_vtbl.h"
#include "RTOp_obj_value_value_vtbl.h"  /* vtbl for operator object instance data */

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )

/* Implementation functions for RTOp_RTOp TOp_force_in_bounds */

static int RTOp_TOp_force_in_bounds_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  RTOp_index_type        sub_dim;
  RTOp_value_type        *x_val;
  ptrdiff_t              x_val_s;
  const RTOp_value_type  *xl_val;
  ptrdiff_t              xl_val_s;
  const RTOp_value_type  *xu_val;
  ptrdiff_t              xu_val_s;
  register RTOp_index_type k;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 2 || vecs == NULL )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 1 || targ_vecs == NULL )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  if( targ_vecs[0].sub_dim != vecs[0].sub_dim
    || targ_vecs[0].sub_dim != vecs[1].sub_dim )
    return RTOp_ERR_INCOMPATIBLE_VECS;

  /* */
  /* Get pointers to data */
  /* */

  /* x */
  sub_dim       = targ_vecs[0].sub_dim;
  x_val         = targ_vecs[0].values;
  x_val_s       = targ_vecs[0].values_stride;
  /* xl */
  xl_val         = vecs[0].values;
  xl_val_s       = vecs[0].values_stride;
  /* xu */
  xu_val         = vecs[1].values;
  xu_val_s       = vecs[1].values_stride;

  /* */
  /* Force in bounds */
  /* */

  if( x_val_s == 1 && xl_val_s == 1 && xu_val_s == 1 ) {
    /* Slightly faster loop for unit stride vectors */
    for( k = 0; k < sub_dim; ++k, ++x_val, ++xl_val, ++xu_val ) {
      if( *x_val < *xl_val )
        *x_val = *xl_val;
      else if( *x_val > *xu_val )
        *x_val = *xu_val;
    }
  }
  else {
    /* More general implementation for non-unit strides */
    for( k = 0; k < sub_dim; ++k, x_val+=x_val_s, xl_val+=xl_val_s, xu_val+=xu_val_s ) {
      if( *x_val < *xl_val )
        *x_val = *xl_val;
      else if( *x_val > *xu_val )
        *x_val = *xu_val;
    }
  }

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_TOp_force_in_bounds_vtbl =
{
  &RTOp_obj_null_vtbl  /* Use a null object for instance data */
  ,&RTOp_obj_null_vtbl /* use null type for target object */
  ,"TOp_force_in_bounds"
  ,NULL /* use default from reduct_vtbl */
  ,RTOp_TOp_force_in_bounds_apply_op
  ,NULL
  ,NULL
};

/* Class specific functions */

int RTOp_TOp_force_in_bounds_construct( struct RTOp_RTOp* op )
{
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_TOp_force_in_bounds_vtbl;
  return 0; /* success? */
}

int RTOp_TOp_force_in_bounds_destroy( struct RTOp_RTOp* op )
{
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0; /* success? */
}

/* Implementation functions for RTOp_RTOp for TOp_force_in_bounds_buffer */

static int RTOp_TOp_force_in_bounds_buffer_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget reduct_obj )
{
  /* */
  /* Declare local variables */
  /* */

  /* Access to the operator object instance data */
  struct RTOp_value_value_type *data_cntr = (struct RTOp_value_value_type*)obj_data;
  RTOp_value_type *rel_push = &data_cntr->value1;
  RTOp_value_type *abs_push = &data_cntr->value2;
  /* Vector data */
  RTOp_index_type           sub_dim;
  /* z0 */
  RTOp_value_type           *z0_val;
  ptrdiff_t                 z0_val_s;
  /* v0 */
  const RTOp_value_type     *v0_val;
  ptrdiff_t                 v0_val_s;
  /* v1 */
  const RTOp_value_type     *v1_val;
  ptrdiff_t                 v1_val_s;

  /* Automatic temporary variables */
  register RTOp_index_type  k;
  /* User defined temporary variables */
  RTOp_value_type xl_sb;
  RTOp_value_type xu_sb;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 2 || ( num_vecs && vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 1 || ( num_targ_vecs && targ_vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  if( /* Validate sub_dim */
    vecs[1].sub_dim != vecs[0].sub_dim
    || targ_vecs[0].sub_dim != vecs[0].sub_dim
    )
    return RTOp_ERR_INCOMPATIBLE_VECS;
  assert(obj_data);


  /* */
  /* Get pointers to data */
  /* */
  sub_dim       = vecs[0].sub_dim;
  /* z0 */
  z0_val        = targ_vecs[0].values;
  z0_val_s      = targ_vecs[0].values_stride;
  /* v0 */
  v0_val        = vecs[0].values;
  v0_val_s      = vecs[0].values_stride;
  /* v1 */
  v1_val        = vecs[1].values;
  v1_val_s      = vecs[1].values_stride;


  /* */
  /* Apply the operator: */
  /* */
  for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s, v1_val += v1_val_s, z0_val += z0_val_s )
  {
    /* Element-wise transformation */
    xl_sb = min((*v0_val) + (*rel_push)*((*v1_val)-(*v0_val)), (*v0_val) + (*abs_push));
    xu_sb = max((*v1_val) - (*rel_push)*((*v1_val)-(*v0_val)), (*v1_val) - (*abs_push));
    if (xl_sb >= xu_sb)
        {
        (*z0_val) = (*v0_val) + ((*v1_val)-(*v0_val))/2.0;
        }
    else if ((*z0_val) < xl_sb)
        { (*z0_val) = xl_sb; }
    else if ((*z0_val) > xu_sb)
        { (*z0_val) = xu_sb; }
    /* Otherwise, leave it */
  }

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_TOp_force_in_bounds_buffer_vtbl =
{
  &RTOp_obj_value_value_vtbl
  ,&RTOp_obj_null_vtbl
  ,"TOp_force_in_bounds_buffer"
  ,NULL
  ,RTOp_TOp_force_in_bounds_buffer_apply_op
  ,NULL
  ,NULL
};

/* Class specific functions */

int RTOp_TOp_force_in_bounds_buffer_construct( RTOp_value_type rel_push, RTOp_value_type abs_push,  struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert(op);
#endif
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_TOp_force_in_bounds_buffer_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  return RTOp_TOp_force_in_bounds_buffer_init(rel_push,abs_push,op);
}

int RTOp_TOp_force_in_bounds_buffer_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0;
}

int RTOp_TOp_force_in_bounds_buffer_init( RTOp_value_type rel_push, RTOp_value_type abs_push, struct RTOp_RTOp* op )
{
  struct RTOp_value_value_type *ptr_data_cntr = (struct RTOp_value_value_type*)op->obj_data;
  RTOp_value_type *ptr_rel_push = &ptr_data_cntr->value1;
  RTOp_value_type *ptr_abs_push = &ptr_data_cntr->value2;
  *ptr_rel_push = rel_push;
  *ptr_abs_push = abs_push;
  return 0;
}

