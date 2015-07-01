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

#include "RTOp_ROp_get_sub_vector.h"
#include "RTOp_obj_free_free.h"

#include <stdlib.h>

#define MY_MIN(a,b) a < b ? a : b

/* Operator object data virtual function table */

struct RTOp_ROp_get_sub_vector_rng_t { /* operator object instance data */
  RTOp_index_type  l;
  RTOp_index_type  u;
};

static int get_op_type_num_entries(
  const struct RTOp_obj_type_vtbl_t*   vtbl
  ,const void* obj_data
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
  assert( num_values );
  assert( num_indexes );
  assert( num_chars );
  *num_values  = 0;
  *num_indexes = 2; /* l, u */
  *num_chars   = 0;
  return 0;
}

static int obj_create(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data
  , RTOp_ReductTarget* obj )
{
	*obj = malloc(sizeof(struct RTOp_ROp_get_sub_vector_rng_t));
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
  const struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  assert(obj_data);
  assert( num_values  == 0 );
  assert( num_indexes == 2 );
  assert( num_chars   == 0 );
    rng = (const struct RTOp_ROp_get_sub_vector_rng_t*)obj_data;
  index_data[0] = rng->l;
  index_data[1] = rng->u;
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
  struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  assert(obj_data);
  assert( num_values  == 0 );
  assert( num_indexes == 2 );
  assert( num_chars   == 0 );
  if(*obj_data == NULL)
    *obj_data = malloc(sizeof(struct RTOp_ROp_get_sub_vector_rng_t));
  rng = (struct RTOp_ROp_get_sub_vector_rng_t*)*obj_data;
  rng->l = index_data[0];
  rng->u = index_data[1];
  return 0;
}

static struct RTOp_obj_type_vtbl_t  instance_obj_vtbl =
{
  get_op_type_num_entries
  ,obj_create
  ,NULL
  ,RTOp_obj_free_free
  ,extract_op_state
  ,load_op_state
};

/* Reduction object virtual function table */

static int get_targ_type_num_entries(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void* obj_data
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
  const struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  assert(obj_data);
  assert( num_values );
  assert( num_indexes );
  assert( num_chars );
    rng = (const struct RTOp_ROp_get_sub_vector_rng_t*)obj_data;
  *num_values  = rng->u - rng->l + 1; /* dense storage for elements of sub-vector to get */
  *num_indexes = 2;                   /* l and u */
  *num_chars   = 0;
  return 0;
}

static int targ_obj_create(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* obj_data
  , RTOp_ReductTarget* targ_obj )
{
  const struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  const int mem_size = sizeof(struct RTOp_SubVector);
  struct RTOp_SubVector *sub_vec_targ = NULL;
  RTOp_index_type sub_dim = 0;
  /* Get the range of the sub-vector */
  assert(obj_data);
  rng = (const struct RTOp_ROp_get_sub_vector_rng_t*)obj_data;
  sub_dim = rng->u - rng->l + 1;
  /* Allocate the sub-vector target object */
  *targ_obj = malloc(mem_size);
  sub_vec_targ = (struct RTOp_SubVector*)*targ_obj;
  /* Setup storage for the target sub-vector */
  RTOp_sub_vector(
     rng->l - 1                                                        /* global_offset */
     ,sub_dim                                                          /* sub_dim */
     ,(const RTOp_value_type*)malloc(sub_dim*sizeof(RTOp_value_type))  /* values[] */
     ,1                                                                /* value_stride */
     ,sub_vec_targ );
  /* Initialize the sub-vector to zero */
  vtbl->obj_reinit( vtbl, obj_data, *targ_obj );
  return 0;
}

static int targ_obj_reinit(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* obj_data
  , RTOp_ReductTarget targ_obj )
{
  const struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  struct RTOp_SubVector *sub_vec_targ = NULL;
  RTOp_index_type sub_dim = 0;
  RTOp_value_type *values = NULL;
  register RTOp_index_type k;
  assert(obj_data);
  /* Get the range of the sub-vector */
    rng = (const struct RTOp_ROp_get_sub_vector_rng_t*)obj_data;
  sub_dim = rng->u - rng->l + 1;
  /* Get the target sub-vector */
  sub_vec_targ = (struct RTOp_SubVector*)targ_obj;
  assert( sub_dim == sub_vec_targ->sub_dim );
  assert( sub_vec_targ->values );
  /* Initialize the values to zero */
  values = (RTOp_value_type*)sub_vec_targ->values;
  for( k = 0; k < sub_dim; ++k )
    *values++ = 0.0;
  return 0;
}

static int targ_obj_free(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* obj_data
  , RTOp_ReductTarget* targ_obj )
{
  const struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  struct RTOp_SubVector *sub_vec_targ = NULL;
  RTOp_index_type sub_dim = 0;
  assert(obj_data);
  /* Get the range of the sub-vector */
    rng = (const struct RTOp_ROp_get_sub_vector_rng_t*)obj_data;
  sub_dim = rng->u - rng->l + 1;
  /* Get the target sub-vector */
  sub_vec_targ = (struct RTOp_SubVector*)*targ_obj;
  assert( sub_dim == sub_vec_targ->sub_dim );
  assert( sub_vec_targ->values );
  /* Deallocate the vectors and the object */
  if( (void*)sub_vec_targ->values )
    free( (void*)sub_vec_targ->values );
  free( (void*)sub_vec_targ );
  *targ_obj = RTOp_REDUCT_OBJ_NULL;
  return 0;
}

static int targ_extract_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void *       obj_data
  ,void *             targ_obj
  ,int                num_values
  ,RTOp_value_type    value_data[]
  ,int                num_indexes
  ,RTOp_index_type    index_data[]
  ,int                num_chars
  ,RTOp_char_type     char_data[]
  )
{
  const struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  struct RTOp_SubVector *sub_vec_targ = NULL;
  RTOp_index_type sub_dim = 0;
  register RTOp_index_type k;
  assert(obj_data);
  /* Get the range of the sub-vector */
  rng = (const struct RTOp_ROp_get_sub_vector_rng_t*)obj_data;
  sub_dim = rng->u - rng->l + 1;
  /* Get the target sub-vector */
  assert( targ_obj  );
  sub_vec_targ = (struct RTOp_SubVector*)targ_obj;
  assert( sub_dim == sub_vec_targ->sub_dim );
  assert( sub_vec_targ->values );
  /* Extract the state */
  assert( num_values  == sub_dim );
  assert( num_indexes == 2 );
  assert( num_chars   == 0 );
  for( k = 0; k < sub_dim; ++k )
    value_data[k] = sub_vec_targ->values[k];
  index_data[0] = rng->l;
  index_data[1] = rng->u;
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
  ,void **                 targ_obj
  )
{
  const struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  struct RTOp_SubVector *sub_vec_targ = NULL;
  RTOp_index_type sub_dim = 0;
  RTOp_value_type *values = NULL;
  register RTOp_index_type k;
  assert(obj_data);
  /* Get the range of the sub-vector */
    rng = (const struct RTOp_ROp_get_sub_vector_rng_t*)obj_data;
  sub_dim = rng->u - rng->l + 1;
  /* Get the target sub-vector */
  assert( *targ_obj );
  sub_vec_targ = (struct RTOp_SubVector*)*targ_obj;
  assert( sub_dim == sub_vec_targ->sub_dim );
  assert( sub_vec_targ->values );
  /* Load the state */
  assert( num_values  == sub_dim );
  assert( num_indexes == 2 );
  assert( num_chars   == 0 );
  assert( index_data[0] == sub_vec_targ->global_offset + 1 );
  assert( index_data[1] == sub_vec_targ->global_offset + sub_vec_targ->sub_dim );
  values = (RTOp_value_type*)sub_vec_targ->values;
  for( k = 0; k < sub_dim; ++k )
    *values++ = value_data[k];
  RTOp_sub_vector(
    rng->l - 1             /* global_offset */
    ,sub_dim               /* sub_dim */
    ,sub_vec_targ->values  /* values[] */
    ,1                     /* value_stide */
    ,sub_vec_targ );
  return 0;
}

static struct RTOp_obj_type_vtbl_t  targ_obj_vtbl =
{
  get_targ_type_num_entries
  ,targ_obj_create
  ,targ_obj_reinit
  ,targ_obj_free
  ,targ_extract_state
  ,targ_load_state
};

/* Implementation functions */

static int RTOp_ROp_get_sub_vector_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  const struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  struct RTOp_SubVector *sub_vec_targ = NULL;
  RTOp_index_type        global_offset;
  RTOp_index_type        sub_dim;
  const RTOp_value_type  *v0_val;
  ptrdiff_t              v0_val_s;
  RTOp_index_type i, i_l, i_u;

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

  /* Get the range of the sub-vector that we are trying to extract */
  assert(obj_data);
  rng = (const struct RTOp_ROp_get_sub_vector_rng_t*)obj_data;

  /* Get the sub-vector target object */
  assert( targ_obj );
  sub_vec_targ = (struct RTOp_SubVector*)targ_obj;
  assert( sub_vec_targ->global_offset + 1 == rng->l );
  assert( sub_vec_targ->global_offset + sub_vec_targ->sub_dim == rng->u );

  /* v0 */
  global_offset  = vecs[0].global_offset;
  sub_dim        = vecs[0].sub_dim;
  v0_val         = vecs[0].values;
  v0_val_s       = vecs[0].values_stride;

  /* */
  /* Extract part of the sub-vector for this chunk. */
  /* */
  /* We only want the elements from rng->l to rng->u */
  /* and this vector chunk only has elements from global_offset */
  /* to global_offset + sub_dim. */
  /* */

  if( rng->u < global_offset + 1 || global_offset + sub_dim < rng->l )
    return 0; /* None of the sub-vector that we are looking for is not in this vector chunk! */

  i_l = ( rng->l <= ( global_offset + 1 )       ? 1       : rng->l - global_offset  );
  i_u = ( rng->u >= ( global_offset + sub_dim ) ? sub_dim : rng->u - global_offset  );

  for( i = i_l; i <= i_u; ++i )
    ((RTOp_value_type*)sub_vec_targ->values)[i-1+(global_offset-(rng->l-1))] = v0_val[(i-1)*v0_val_s];

  return 0; /* success? */

}

static int reduce_reduct_objs(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data /* Can be NULL! */
  , RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj )
{
  const struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  struct RTOp_SubVector
    *in_sub_vec_targ    = NULL,
    *inout_sub_vec_targ = NULL;
  RTOp_index_type sub_dim = 0;
  RTOp_value_type *inout_values = NULL;
  register RTOp_index_type k;
  /* Get the range of the sub-vector */
  assert(obj_data);
    rng = (const struct RTOp_ROp_get_sub_vector_rng_t*)obj_data;
  sub_dim = rng->u - rng->l + 1;
  /* Get the in target sub-vector */
  assert( in_reduct_obj );
  in_sub_vec_targ = (struct RTOp_SubVector*)in_reduct_obj;
  assert( sub_dim == in_sub_vec_targ->sub_dim );
  assert( in_sub_vec_targ->values );
  /* Get the inout target sub-vector */
  assert( inout_reduct_obj );
  inout_sub_vec_targ = (struct RTOp_SubVector*)inout_reduct_obj;
  assert( sub_dim == inout_sub_vec_targ->sub_dim );
  assert( inout_sub_vec_targ->values );
  /* Perform the reduction */
  inout_values = (RTOp_value_type*)inout_sub_vec_targ->values;
  for( k = 0; k < sub_dim; ++k )
    *inout_values++ += in_sub_vec_targ->values[k];
  return 0;
}

static void CALL_API external_reduct_op( void* in_targ_array, void* inout_targ_array
  , int* len, RTOp_Datatype* datatype )
{
  int num_values_off, num_indexes_off, num_chars_off
    , values_off, l_off, u_off;
  int num_values, num_indexes, num_chars;
  const RTOp_value_type *in_values    = NULL;
  RTOp_value_type       *inout_values = NULL;
  register RTOp_index_type k;
  /* Get the offsets for the number of elements of each data type */
  num_values_off  = 0;
  num_indexes_off = num_values_off  + sizeof(RTOp_value_type);
  num_chars_off   = num_indexes_off + sizeof(RTOp_value_type);
  /* Get the number of elements of each data type */
  num_values  = *(RTOp_value_type*)((char*)in_targ_array + num_values_off);
  num_indexes = *(RTOp_value_type*)((char*)in_targ_array + num_indexes_off);
  num_chars   = *(RTOp_value_type*)((char*)in_targ_array + num_chars_off);
#ifdef RTOp_DEBUG
  assert( num_values > 0 );
  assert( num_indexes == 2 );
  assert( num_chars   == 0 );
  assert( num_values  == *(RTOp_value_type*)((char*)inout_targ_array + num_values_off) );
  assert( num_indexes == *(RTOp_value_type*)((char*)inout_targ_array + num_indexes_off) );
  assert( num_chars   == *(RTOp_value_type*)((char*)inout_targ_array + num_chars_off) );
#endif
  /* Get the offsets for the sub-vector values and range l and u */
  values_off = num_chars_off + sizeof(RTOp_value_type);
  l_off      = values_off    + num_values * sizeof(RTOp_value_type);
  u_off      = l_off         + sizeof(RTOp_index_type);
#ifdef RTOp_DEBUG
  assert( *(RTOp_index_type*)((char*)in_targ_array + l_off)
      == *(RTOp_index_type*)((char*)inout_targ_array + l_off) );
  assert( *(RTOp_index_type*)((char*)in_targ_array + u_off)
      == *(RTOp_index_type*)((char*)inout_targ_array + u_off) );
#endif
  /* Perform the reduction! (just add the elements together) */
  in_values    = (const RTOp_value_type*)((char*)in_targ_array + values_off);
  inout_values = (RTOp_value_type*)((char*)inout_targ_array + values_off);
  for( k = 0; k < num_values; ++k )
    *inout_values++ += *in_values++;
}

static int get_reduct_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr )
{
  *reduct_op_func_ptr = external_reduct_op;
  return 0;
}

/* Public interface */

const struct RTOp_RTOp_vtbl_t RTOp_ROp_get_sub_vector_vtbl =
{
  &instance_obj_vtbl
  ,&targ_obj_vtbl
  ,"RTOp_ROp_get_sub_vector"
  ,NULL
  ,RTOp_ROp_get_sub_vector_apply_op
  ,reduce_reduct_objs
  ,get_reduct_op
};

int RTOp_ROp_get_sub_vector_construct(
  RTOp_index_type l, RTOp_index_type u, struct RTOp_RTOp* op )
{
  struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  op->vtbl = &RTOp_ROp_get_sub_vector_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  rng = (struct RTOp_ROp_get_sub_vector_rng_t*)op->obj_data;
  rng->l = l;
  rng->u = u;
  return 0;
}

int RTOp_ROp_get_sub_vector_set_range(
  RTOp_index_type l, RTOp_index_type u, struct RTOp_RTOp* op )
{
  struct RTOp_ROp_get_sub_vector_rng_t *rng = NULL;
  assert( op->vtbl );
  assert( op->obj_data );
  rng = (struct RTOp_ROp_get_sub_vector_rng_t*)op->obj_data;
  rng->l = l;
  rng->u = u;
  return 0;
}

int RTOp_ROp_get_sub_vector_destroy( struct RTOp_RTOp* op )
{
#ifdef TEUCHOS_DEBUG
  assert( op->vtbl );
  assert( op->obj_data );
#endif
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->vtbl = NULL;
  return 0;
}

struct RTOp_SubVector RTOp_ROp_get_sub_vector_val(
  RTOp_ReductTarget targ_obj
  )
{
  return *((struct RTOp_SubVector*)targ_obj);
}
