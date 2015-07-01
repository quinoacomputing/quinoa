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

#ifndef RTOP_GET_REDUCT_OP_H
#define RTOP_GET_REDUCT_OP_H

#define INSERT_GET_REDUCT_OP_FUNCS(                                                                                   \
    num_values,num_indexes,num_chars,reduct_obj_t,reduce_reduct_obj,targ_load_state,targ_extract_state                \
    ,external_reduct_op,get_reduct_op                                                                                 \
    )                                                                                                                 \
static void CALL_API external_reduct_op( void* in_targ_array, void* inout_targ_array                                           \
  , int* len, RTOp_Datatype* datatype )                                                                             \
{                                                                                                                     \
  struct reduct_obj_t                                                                                               \
    in_obj, *in_obj_p = &in_obj, inout_obj, *inout_obj_p = &inout_obj;                                            \
    const int                                                                                                         \
        values_off  = 3*sizeof(RTOp_value_type),                                                                      \
        indexes_off = values_off + num_values*sizeof(RTOp_value_type),                                                \
        chars_off   = indexes_off + num_indexes*sizeof(RTOp_index_type);                                              \
  const int size_obj = chars_off + (num_chars)*sizeof(RTOp_index_type);                                             \
  char                                                                                                              \
     *in_array    = in_targ_array,                                                                                 \
    *inout_array = inout_targ_array;                                                                              \
  int i;                                                                                                            \
  for( i = 0; i < *len; ++i, in_array += size_obj, inout_array += size_obj ) {                                      \
    targ_load_state(                                                                                              \
      NULL, NULL                                                                                                \
      ,num_values,  num_values  ? (RTOp_value_type*)(in_array + values_off) : NULL                              \
      ,num_indexes, num_indexes ? (RTOp_index_type*)(in_array + values_off) : NULL                              \
      ,num_chars,   num_chars   ? (RTOp_char_type*) (in_array + values_off) : NULL                              \
      ,(void**)&in_obj_p );                                                                                     \
    targ_load_state(                                                                                              \
      NULL, NULL                                                                                                \
      ,num_values,  num_values  ? (RTOp_value_type*)(inout_array + values_off) : NULL                           \
      ,num_indexes, num_indexes ? (RTOp_index_type*)(inout_array + values_off) : NULL                           \
      ,num_chars,   num_chars   ? (RTOp_char_type*) (inout_array + values_off) : NULL                           \
      ,(void**)&inout_obj_p );                                                                                  \
    reduce_reduct_objs( NULL, NULL, &in_obj, &inout_obj );                                                        \
    targ_extract_state(                                                                                           \
      NULL, NULL, &inout_obj                                                                                    \
      ,num_values,  num_values  ? (RTOp_value_type*)(inout_array + values_off) : NULL                           \
      ,num_indexes, num_indexes ? (RTOp_index_type*)(inout_array + values_off) : NULL                           \
      ,num_chars,   num_chars   ? (RTOp_char_type*) (inout_array + values_off) : NULL                           \
      );                                                                                                        \
  }                                                                                                                 \
}                                                                                                                     \
static int get_reduct_op(                                                                                             \
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data                                                         \
  , RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr )                                                                 \
{                                                                                                                     \
  *reduct_op_func_ptr = external_reduct_op;                                                                         \
  return 0;                                                                                                         \
}

#endif
