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

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "RTOp.h"

/* Misc. */

void RTOp_sub_vector(
  RTOp_index_type global_offset, RTOp_index_type sub_dim
  ,const RTOp_value_type values[], ptrdiff_t values_stride
  ,struct RTOp_SubVector *sub_vec
  )
{
  /* Validate input */
#ifdef RTOp_DEBUG
  assert( values != NULL );
#endif
  /* Set members */
  sub_vec->global_offset  = global_offset;
  sub_vec->sub_dim        = sub_dim;
  sub_vec->values         = values;
  sub_vec->values_stride  = values_stride;
}

void RTOp_sub_vector_null( struct RTOp_SubVector *sub_vec )
{
  sub_vec->global_offset  = 0;
  sub_vec->sub_dim        = 0;
  sub_vec->values         = NULL;
  sub_vec->values_stride  = 0;
}

void RTOp_mutable_sub_vector(
  RTOp_index_type global_offset, RTOp_index_type sub_dim
  ,RTOp_value_type values[], ptrdiff_t values_stride
  ,struct RTOp_MutableSubVector *sub_vec
  )
{
  /* Validate input */
#ifdef RTOp_DEBUG
  assert( sub_vec );
  assert( values != NULL );
#endif
  /* Set members */
  sub_vec->global_offset  = global_offset;
  sub_vec->sub_dim        = sub_dim;
  sub_vec->values         = values;
  sub_vec->values_stride  = values_stride;
}

void RTOp_mutable_sub_vector_null( struct RTOp_MutableSubVector *sub_vec )
{
#ifdef RTOp_DEBUG
  assert( sub_vec );
#endif
  sub_vec->global_offset  = 0;
  sub_vec->sub_dim        = 0;
  sub_vec->values         = NULL;
  sub_vec->values_stride  = 0;
}

/* RTOp_RTOp */

int RTOp_get_op_name(
	const struct RTOp_RTOp* op
	,const char** op_name
	)
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->op_name );
  assert( op_name );
#endif
  *op_name = op->vtbl->op_name;
  return 0;
}

int RTOp_get_op_type_num_entries(
	const struct RTOp_RTOp* op
	,int* num_values
	,int* num_indexes
	,int* num_chars
	)
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->obj_data_vtbl );
  assert( op->vtbl->obj_data_vtbl->get_obj_type_num_entries );
  assert( num_values );
  assert( num_indexes );
  assert( num_chars );
#endif
  return op->vtbl->obj_data_vtbl->get_obj_type_num_entries(
    op->vtbl->obj_data_vtbl,op->obj_data
    , num_values, num_indexes, num_chars );
}

int RTOp_extract_op_state(
  const struct RTOp_RTOp* op
  ,int num_values
  ,RTOp_value_type value_data[]
  ,int num_indexes
  ,RTOp_index_type index_data[]
  ,int num_chars
  ,RTOp_char_type  char_data[]
  )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->obj_data_vtbl );
  assert( op->vtbl->obj_data_vtbl->extract_state );
  if(num_values)  assert( value_data );
  if(num_indexes) assert( index_data );
  if(num_chars)   assert( char_data );
#endif
  return op->vtbl->obj_data_vtbl->extract_state(
    op->vtbl->obj_data_vtbl, NULL, op->obj_data
    , num_values, value_data
    , num_indexes, index_data
    , num_chars, char_data );
}

int RTOp_load_op_state(
  int num_values
  ,const RTOp_value_type value_data[]
  ,int num_indexes
  ,const RTOp_index_type index_data[]
  ,int num_chars
  ,const RTOp_char_type  char_data[]
  ,struct RTOp_RTOp* op
  )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->obj_data_vtbl );
  assert( op->vtbl->obj_data_vtbl->load_state );
  if(num_values)  assert( value_data );
  if(num_indexes) assert( index_data );
  if(num_chars)   assert( char_data );
#endif
  return op->vtbl->obj_data_vtbl->load_state(
    op->vtbl->obj_data_vtbl, NULL
    ,num_values,  value_data
    ,num_indexes, index_data
    ,num_chars,   char_data
    , &op->obj_data );
}

int RTOp_free_op( struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->obj_data_vtbl );
  assert( op->vtbl->obj_data_vtbl->obj_free );
#endif
  return op->vtbl->obj_data_vtbl->obj_free(
    op->vtbl->obj_data_vtbl,NULL,&op->obj_data );
}

int RTOp_get_reduct_type_num_entries(
  const struct RTOp_RTOp* op
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->reduct_vtbl );
  assert( op->vtbl->reduct_vtbl->get_obj_type_num_entries );
  assert( num_values );
  assert( num_indexes );
  assert( num_chars );
#endif
  return op->vtbl->reduct_vtbl->get_obj_type_num_entries(
    op->vtbl->reduct_vtbl, op->obj_data
    , num_values, num_indexes, num_chars );
}

int RTOp_reduct_obj_create( const struct RTOp_RTOp* op
  , RTOp_ReductTarget* reduct_obj )
{
  int err = 0;
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->reduct_vtbl );
  assert( op->vtbl->reduct_vtbl->obj_create );
#endif
  err = op->vtbl->reduct_vtbl->obj_create(
    op->vtbl->reduct_vtbl,op->obj_data,reduct_obj );
  if(err) return err;
  if( op->vtbl->reduct_obj_reinit )
    return op->vtbl->reduct_obj_reinit(
      op->vtbl,op->obj_data,*reduct_obj );
  return err;
}

int RTOp_reduct_obj_reinit( const struct RTOp_RTOp* op
  , RTOp_ReductTarget reduct_obj )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->reduct_vtbl );
  assert( op->vtbl->reduct_vtbl->obj_reinit );
  assert( reduct_obj != RTOp_REDUCT_OBJ_NULL );
#endif
  if( op->vtbl->reduct_obj_reinit )
    return op->vtbl->reduct_obj_reinit(
      op->vtbl,op->obj_data,reduct_obj );
  else
    return op->vtbl->reduct_vtbl->obj_reinit(
      op->vtbl->reduct_vtbl,op->obj_data,reduct_obj );
}

int RTOp_reduct_obj_free( const struct RTOp_RTOp* op
  , RTOp_ReductTarget* reduct_obj )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->reduct_vtbl );
  assert( op->vtbl->reduct_vtbl->obj_free );
#endif
  return op->vtbl->reduct_vtbl->obj_free(
    op->vtbl->reduct_vtbl,op->obj_data,reduct_obj);
}

int RTOp_extract_reduct_obj_state(
  const struct RTOp_RTOp*   op
  ,const RTOp_ReductTarget  reduct_obj
  ,int                      num_values
  ,RTOp_value_type          value_data[]
  ,int                      num_indexes
  ,RTOp_index_type          index_data[]
  ,int                      num_chars
  ,RTOp_char_type           char_data[]
  )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->reduct_vtbl );
  assert( op->vtbl->reduct_vtbl->extract_state );
  if(num_values)  assert( value_data );
  if(num_indexes) assert( index_data );
  if(num_chars)   assert( char_data );
#endif
  return op->vtbl->reduct_vtbl->extract_state(
    op->vtbl->reduct_vtbl, op->obj_data, reduct_obj
    ,num_values,  value_data
    ,num_indexes, index_data
    ,num_chars,   char_data );
}

int RTOp_load_reduct_obj_state(
  const struct RTOp_RTOp*  op
  ,int                     num_values
  ,const RTOp_value_type   value_data[]
  ,int                     num_indexes
  ,const RTOp_index_type   index_data[]
  ,int                     num_chars
  ,const RTOp_char_type    char_data[]
  ,RTOp_ReductTarget       reduct_obj
  )
{
  RTOp_ReductTarget reduct_obj_in = reduct_obj; /* Just keep a reference */
  int err = 0;
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->reduct_vtbl );
  assert( op->vtbl->reduct_vtbl->load_state );
  if(num_values)  assert( value_data );
  if(num_indexes) assert( index_data );
  if(num_chars)   assert( char_data );
#endif
  err = op->vtbl->reduct_vtbl->load_state(
    op->vtbl->reduct_vtbl, op->obj_data
    ,num_values,  value_data
    ,num_indexes, index_data
    ,num_chars,   char_data
    ,&reduct_obj );
  if(err != 0) return err;
  if( reduct_obj != reduct_obj_in )
    return RTOp_ERR_INVALID_USAGE;
  return 0; /* Success! */
}

int RTOp_coord_invariant(
	const struct RTOp_RTOp   *op
	,int                     *coord_invariant
	)
{
#ifdef RTOp_DEBUG
	assert( op );
	assert( op->vtbl );
/*	assert( op->vtbl->coord_invariant ); */
	assert( coord_invariant );
#endif
/*	return op->vtbl->coord_invariant(op->vtbl,op->obj_data,coord_invariant); */
	*coord_invariant = 1; /* ToDo: Implement the above code! */
	return 0;
}

int RTOp_apply_op( const struct RTOp_RTOp* op
  , const int num_vecs, const struct RTOp_SubVector sub_vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_sub_vecs[]
  , RTOp_ReductTarget reduct_obj )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->apply_op );
  if(num_vecs) assert(sub_vecs);
  if(num_targ_vecs) assert(targ_sub_vecs);
#endif
  return op->vtbl->apply_op(
    op->vtbl,op->obj_data,num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs
    ,reduct_obj);
}

int RTOp_reduce_reduct_objs( const struct RTOp_RTOp* op
  , RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->reduce_reduct_objs );
#endif
  return op->vtbl->reduce_reduct_objs(op->vtbl,op->obj_data,in_reduct_obj,inout_reduct_obj);
}

int RTOp_get_reduct_op( const struct RTOp_RTOp* op
  , RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr )
{
#ifdef RTOp_DEBUG
  assert( op );
  assert( op->vtbl );
  assert( op->vtbl->get_reduct_op );
#endif
  return op->vtbl->get_reduct_op(op->vtbl,op->obj_data,reduct_op_func_ptr);
}

/* */
/* RTOp_Server */
/* */

/* RTOp_Server data */

static int RTOp_Server_num_ops = 0;

#define RTOp_SERVER_MAX_NUM_ENTRIES 50
struct RTOp_Server_op_class_name {
  char name[RTOp_SERVER_MAX_NUM_ENTRIES+1];
};
static struct RTOp_Server_op_class_name
  RTOp_Server_op_names[RTOp_SERVER_MAX_NUM_ENTRIES];
typedef const struct RTOp_RTOp_vtbl_t* RTOp_RTOp_vtbl_t_ptr;
static RTOp_RTOp_vtbl_t_ptr
  RTOp_Server_op_vtbl[RTOp_SERVER_MAX_NUM_ENTRIES];

/* */
/* Private RTOp_Server functions */
/* */
/* In this simplistic implementation the names and vtbl pointers */
/* are stored an accessed in unsorted tables. */
/* */

/* returns the position in the table where this name is found. */
/* otherwise it returns < 0. */
static int find_op_name( const struct RTOp_Server_op_class_name op_name_tbl[]
  , int num_entries, const char op_name[] )
{
  int k = 0;
  for( k = 0; k < num_entries; ++k )
  {
    if( strcmp( op_name_tbl[k].name, op_name ) == 0 )
      return k;
  }
  return -1; /* Did not find the name */
}

/* returns the position in the table where this reduct vtbl is found. */
/* otherwise it returns < 0. */
static int find_op_vtbl( const RTOp_RTOp_vtbl_t_ptr op_vtbl_tbl[]
  , int num_entries, RTOp_RTOp_vtbl_t_ptr op_vtbl )
{
  int k = 0;
  for( k = 0; k < num_entries; ++k )
  {
    if( op_vtbl_tbl[k] == op_vtbl )
      return k;
  }
  return -1; /* Did not find the name */
}

/* */
/* Public RTOp_Server functions */
/* */

int RTOp_Server_add_op_name_vtbl( const char op_class_name[]
  , const struct RTOp_RTOp_vtbl_t* op_class_vtbl )
{
  int k = 0;
  if( strlen( op_class_name ) > RTOp_MAX_REDUCT_TRANS_OP_CLASS_NAME )
    return RTOp_SERVER_OP_NAME_TOO_LONG;
  if( (k = find_op_name( RTOp_Server_op_names
    , RTOp_Server_num_ops, op_class_name ) ) >= 0 )
  {
    if( RTOp_Server_op_vtbl[k] != op_class_vtbl )
      return RTOp_SERVER_INCOMPATIBLE_OPS;
    return k; /* The name exits but vtble is the same */
  }
  strcpy( RTOp_Server_op_names[RTOp_Server_num_ops].name, op_class_name );
  RTOp_Server_op_vtbl[RTOp_Server_num_ops] = op_class_vtbl;
  ++RTOp_Server_num_ops;
  return 0; /* Successfully added it */
}

int RTOp_Server_lookup_op_name( const struct RTOp_RTOp_vtbl_t* op_class_vtbl
  , char op_class_name[] )
{
  int k = 0;
  if( ( k = find_op_vtbl(RTOp_Server_op_vtbl,RTOp_Server_num_ops,op_class_vtbl) ) >= 0 )
  {
    strcpy( op_class_name, RTOp_Server_op_names[k].name );
    return 0; /* Success */
  }
  return -1; /* Did not find vtbl */
}

int RTOp_Server_construct_op(
  const char op_class_name[]
  ,int num_values
  ,const RTOp_value_type value_data[]
  ,int num_indexes
  ,const RTOp_index_type index_data[]
  ,int num_chars
  ,const RTOp_char_type  char_data[]
  ,struct RTOp_RTOp* op
  )
{
  int k = 0;
  int err = 0; /* success? */
  if( strlen( op_class_name ) > RTOp_MAX_REDUCT_TRANS_OP_CLASS_NAME )
    return RTOp_SERVER_OP_NAME_TOO_LONG;
  if( ( k = find_op_name(RTOp_Server_op_names,RTOp_Server_num_ops, op_class_name) ) >= 0 )
  {
    op->obj_data = NULL; /* Will be dyn allocated below! */
    op->vtbl     = RTOp_Server_op_vtbl[k];
    err = RTOp_load_op_state(
      num_values,value_data,num_indexes,index_data,num_chars,char_data
      ,op);
    return err;
  }
  return -1; /* op_class_name not found */
}

void RTOp_Server_dump( FILE* file )
{
  int k = 0;
  int jn = 0, jv = 0;
  fprintf( file, "Class names and vtbl pointers for RTOp_RTOp subcasses\n" );
  for( k = 0; k < RTOp_Server_num_ops; ++k ) {
    jn = find_op_name( RTOp_Server_op_names
    , RTOp_Server_num_ops, RTOp_Server_op_names[k].name );
    jv = find_op_vtbl( RTOp_Server_op_vtbl, RTOp_Server_num_ops
      , RTOp_Server_op_vtbl[k] );
    fprintf( file
      , "    class name            = \"%s\"\n"
        "    looked up class name  = \"%s\"\n"
        "    vtbl                  = %p\n"
        "    looked up vtbl        = %p\n"
      , RTOp_Server_op_names[k].name
      , RTOp_Server_op_names[jn].name
      , RTOp_Server_op_vtbl[k]
      , RTOp_Server_op_vtbl[jv]
      );
  }
}
