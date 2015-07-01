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

#include "RTOp_TOp_set_sub_vector.h"
#include "RTOp_obj_null_vtbl.h"

#include <stdlib.h>

/* Operator object data virtual function table */

struct RTOp_TOp_set_sub_vector_state_t { /* operator object instance data */
  struct RTOp_SparseSubVector  sub_vec;   /* The sub vector to set! */
  int                    owns_mem;  /* if true then memory is sub_vec must be deleted! */
};

static int get_op_type_num_entries(
  const struct RTOp_obj_type_vtbl_t*   vtbl
  ,const void* obj_data
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
  const struct RTOp_TOp_set_sub_vector_state_t *state = NULL;
#ifdef RTOp_DEBUG
  assert( num_values );
  assert( num_indexes );
  assert( num_chars );
  assert( obj_data );
#endif
  state = (const struct RTOp_TOp_set_sub_vector_state_t*)obj_data;
  *num_values  = state->sub_vec.sub_nz;  /* values[] */
  *num_indexes =
    8  /* global_offset, sub_dim, sub_nz, values_stride, indices_stride, local_offset, is_sorted, owns_mem */
    + (state->sub_vec.indices ? state->sub_vec.sub_nz : 0 ); /* indices[] */
  *num_chars   = 0;
  return 0;
}

static int op_create(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data
  , RTOp_ReductTarget* obj )
{
  struct RTOp_TOp_set_sub_vector_state_t *state = NULL;
  *obj = malloc(sizeof(struct RTOp_TOp_set_sub_vector_state_t));
  state = (struct RTOp_TOp_set_sub_vector_state_t*)*obj;
  RTOp_sparse_sub_vector_null( &state->sub_vec );
  state->owns_mem = 0;
  return 0;
}

static int op_free(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* dummy
  , void ** obj_data )
{
  struct RTOp_TOp_set_sub_vector_state_t *state = NULL;
#ifdef RTOp_DEBUG
  assert( *obj_data );
#endif
  state = (struct RTOp_TOp_set_sub_vector_state_t*)*obj_data;
  if( state->owns_mem ) {
    free( (void*)state->sub_vec.values );
    free( (void*)state->sub_vec.indices );
  }
  free( *obj_data );
  *obj_data = NULL;
  return 0;
}

static int extract_op_state(
  const struct RTOp_obj_type_vtbl_t*   vtbl
  ,const void *       dummy
  ,void *             obj_data
  ,int                num_values
  ,RTOp_value_type    value_data[]
  ,int                num_indexes
  ,RTOp_index_type    index_data[]
  ,int                num_chars
  ,RTOp_char_type     char_data[]
  )
{
  const struct RTOp_TOp_set_sub_vector_state_t *state = NULL;
  register RTOp_index_type k, j;
#ifdef RTOp_DEBUG
  assert( obj_data );
#endif
  state = (const struct RTOp_TOp_set_sub_vector_state_t*)obj_data;
#ifdef RTOp_DEBUG
  assert( num_values == state->sub_vec.sub_nz );
  assert( num_indexes == 8 + (state->sub_vec.indices ? state->sub_vec.sub_nz : 0 ) );
  assert( num_chars == 0 );
  assert( value_data );
  assert( index_data );
#endif
  for( k = 0; k < state->sub_vec.sub_nz; ++k )
    value_data[k] = state->sub_vec.values[k];
  index_data[k=0] = state->sub_vec.global_offset;
  index_data[++k] = state->sub_vec.sub_dim;
  index_data[++k] = state->sub_vec.sub_nz;
  index_data[++k] = state->sub_vec.values_stride;
  index_data[++k] = state->sub_vec.indices_stride;
  index_data[++k] = state->sub_vec.local_offset;
  index_data[++k] = state->sub_vec.is_sorted;
  index_data[++k] = state->owns_mem;
  if( state->sub_vec.indices ) {
    for( j = 0; j < state->sub_vec.sub_nz; ++j )
      index_data[++k] = state->sub_vec.indices[j];
  }
  return 0;
}

static int load_op_state(
  const struct RTOp_obj_type_vtbl_t*   vtbl
  ,const void *            dummy
  ,int                     num_values
  ,const RTOp_value_type   value_data[]
  ,int                     num_indexes
  ,const RTOp_index_type   index_data[]
  ,int                     num_chars
  ,const RTOp_char_type    char_data[]
  ,void **                 obj_data
  )
{
  struct RTOp_TOp_set_sub_vector_state_t *state = NULL;
  register RTOp_index_type k, j;
  RTOp_value_type *values = NULL;
  /* Allocate the operator object's state data if it has not been already */
  if( *obj_data == NULL ) {
    *obj_data = (void*)malloc(sizeof(struct RTOp_TOp_set_sub_vector_state_t));
    state = (struct RTOp_TOp_set_sub_vector_state_t*)*obj_data;
    RTOp_sparse_sub_vector_null( &state->sub_vec );
    state->owns_mem = 0;
  }
  /* Get the operator object's state data */
#ifdef RTOp_DEBUG
  assert( *obj_data );
#endif
  state = (struct RTOp_TOp_set_sub_vector_state_t*)*obj_data;
#ifdef RTOp_DEBUG
  if( num_indexes > 8 ) assert( num_values == num_indexes - 8 );
#endif
  /* Reallocate storage if we have to */
  if( num_values != state->sub_vec.sub_nz || !state->owns_mem ) {
    /* Delete current storage if owned */
    if( state->owns_mem ) {
      free( (RTOp_value_type*)state->sub_vec.values );
      free( (RTOp_index_type*)state->sub_vec.indices );
    }
    /* We need to reallocate storage for values[] and perhaps */
    state->sub_vec.values = (RTOp_value_type*)malloc(num_values*sizeof(RTOp_value_type));
    if( num_indexes > 8 )
      state->sub_vec.indices = (RTOp_index_type*)malloc(num_values*sizeof(RTOp_index_type));
    state->owns_mem = 1;
  }
  /* Set the internal state! */
#ifdef RTOp_DEBUG
  assert( num_chars == 0 );
  assert( value_data );
  assert( index_data );
#endif
  for( values = (RTOp_value_type*)state->sub_vec.values, k = 0; k < num_values; ++k )
    *values++ = value_data[k];
  state->sub_vec.global_offset  = index_data[k=0];
  state->sub_vec.sub_dim        = index_data[++k];
  state->sub_vec.sub_nz         = index_data[++k];
  state->sub_vec.values_stride  = index_data[++k];
  state->sub_vec.indices_stride = index_data[++k];
  state->sub_vec.local_offset   = index_data[++k];
  state->sub_vec.is_sorted      = index_data[++k];
  state->owns_mem               = index_data[++k];
  if( num_indexes > 8 ) {
    for( j = 0; j < num_values; ++j )
      ((RTOp_index_type*)state->sub_vec.indices)[j] = index_data[++k];
  }
  return 0;
}

static struct RTOp_obj_type_vtbl_t  instance_obj_vtbl =
{
  get_op_type_num_entries
  ,op_create
  ,NULL
  ,op_free
  ,extract_op_state
  ,load_op_state
};

/* Implementation functions */

static int RTOp_TOp_set_sub_vector_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  const struct RTOp_TOp_set_sub_vector_state_t *state = NULL;
  RTOp_index_type        v_global_offset;
  RTOp_index_type        v_sub_dim;
  RTOp_index_type        v_sub_nz;
  const RTOp_value_type  *v_val;
  ptrdiff_t              v_val_s;
  const RTOp_index_type  *v_ind;
  ptrdiff_t              v_ind_s;
  ptrdiff_t              v_l_off;
  int                    v_sorted;
  RTOp_index_type        z_global_offset;
  RTOp_index_type        z_sub_dim;
  RTOp_value_type        *z_val;
  ptrdiff_t              z_val_s;
  register RTOp_index_type  k, i;
  RTOp_index_type num_overlap;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 1 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  assert(targ_vecs);
  assert(targ_obj == RTOp_REDUCT_OBJ_NULL);

  /* */
  /* Get pointers to data */
  /* */

  /* Get the sub-vector we are reading from */
  assert(obj_data);
  state = (const struct RTOp_TOp_set_sub_vector_state_t*)obj_data;

  /* v (the vector to read from) */
  v_global_offset  = state->sub_vec.global_offset;
  v_sub_dim        = state->sub_vec.sub_dim;
  v_sub_nz         = state->sub_vec.sub_nz;
  v_val            = state->sub_vec.values;
  v_val_s          = state->sub_vec.values_stride;
  v_ind            = state->sub_vec.indices;
  v_ind_s          = state->sub_vec.indices_stride;
  v_l_off          = state->sub_vec.local_offset;
  v_sorted         = state->sub_vec.is_sorted;

  /* z (the vector to set) */
  z_global_offset  = targ_vecs[0].global_offset;
  z_sub_dim        = targ_vecs[0].sub_dim;
  z_val            = targ_vecs[0].values;
  z_val_s          = targ_vecs[0].values_stride;

  /* */
  /* Set part of the sub-vector for this chunk. */
  /* */

  if( v_global_offset + v_sub_dim < z_global_offset + 1
    || z_global_offset + z_sub_dim < v_global_offset + 1 )
    return 0; /* The sub-vector that we are setting does not overlap with this vector chunk! */

  if( v_sub_nz == 0 )
    return 0; /* The sparse sub-vector we are reading from is empty? */

  /* Get the number of elements that overlap */
  if( v_global_offset <= z_global_offset ) {
    if( v_global_offset + v_sub_dim >= z_global_offset + z_sub_dim )
      num_overlap = z_sub_dim;
    else
      num_overlap = (v_global_offset + v_sub_dim) - z_global_offset;
  }
  else {
    if( z_global_offset + z_sub_dim >= v_global_offset + v_sub_dim )
      num_overlap = v_sub_dim;
    else
      num_overlap = (z_global_offset + z_sub_dim) - v_global_offset;
  }

  /* Set the part of the sub-vector that overlaps */
  if( v_ind != NULL ) {
    /* Sparse elements */
    /* Set the overlapping elements to zero first. */
    if( v_global_offset >= z_global_offset )
      z_val += (v_global_offset - z_global_offset) * z_val_s;
    for( k = 0; k < num_overlap; ++k, z_val += z_val_s )
      *z_val = 0.0;
    /* Now set the sparse entries */
    z_val = targ_vecs[0].values;
    for( k = 0; k < v_sub_nz; ++k, v_val += v_val_s, v_ind += v_ind_s ) {
      i = v_global_offset + v_l_off + (*v_ind);
      if( z_global_offset < i && i <= z_global_offset + z_sub_dim )
        z_val[ z_val_s * (i - z_global_offset - 1) ] = *v_val;
    }
    /* ToDo: Implement a faster version for v sorted and eliminate the */
    /* if statement in the loop. */
  }
  else {
    /* Dense elemements */
    if( v_global_offset <= z_global_offset )
      v_val += (z_global_offset - v_global_offset) * v_val_s;
    else
      z_val += (v_global_offset - z_global_offset) * z_val_s;
    for( k = 0; k < num_overlap; ++k, v_val += v_val_s, z_val += z_val_s )
      *z_val = *v_val;
  }

  return 0; /* success? */
}

/* Public interface */

const struct RTOp_RTOp_vtbl_t RTOp_TOp_set_sub_vector_vtbl =
{
  &instance_obj_vtbl
  ,&RTOp_obj_null_vtbl /* use null type for reduction target object */
  ,"RTOp_TOp_set_sub_vector"
  ,NULL /* use default from reduct_vtbl */
  ,RTOp_TOp_set_sub_vector_apply_op
  ,NULL /* reduce_reduct_obj */
  ,NULL /* get_reduct_op */
};

int RTOp_TOp_set_sub_vector_construct(
  const struct RTOp_SparseSubVector* sub_vec, struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert(sub_vec);
  assert(op);
#endif
  op->vtbl     = &RTOp_TOp_set_sub_vector_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  return RTOp_TOp_set_sub_vector_set_sub_vec(sub_vec,op);
}

int RTOp_TOp_set_sub_vector_set_sub_vec(
  const struct RTOp_SparseSubVector* sub_vec, struct RTOp_RTOp* op )
{
  struct RTOp_TOp_set_sub_vector_state_t *state = NULL;
#ifdef RTOp_DEBUG
  assert( op->obj_data );
#endif
  state = (struct RTOp_TOp_set_sub_vector_state_t*)op->obj_data;
  if( state->owns_mem ) {
    free( (void*)state->sub_vec.values );
    free( (void*)state->sub_vec.indices );
  }
  state->sub_vec = *sub_vec;  /* We do not own the arrays values[] and indices[] */
  state->owns_mem = 0;
  return 0;
}

int RTOp_TOp_set_sub_vector_destroy( struct RTOp_RTOp* op )
{
  op_free(NULL,NULL,&op->obj_data);
  op->vtbl = NULL;
  return 0;
}
