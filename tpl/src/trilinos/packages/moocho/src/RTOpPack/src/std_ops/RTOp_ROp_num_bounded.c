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

#include "RTOp_ROp_num_bounded.h"
#include "RTOp_obj_value_vtbl.h"
#include "RTOp_reduct_sum_value.h"

/* Note that the reduction quantity that we are accumulating (num_bounded) */
/* is an integral type and really should be delcared as RTOp_index_type. */
/* However, the machinary is already there for accumulating an RTOp_value_type */
/* reduction object so this implementation is just lazy and uses a double */
/* for an integer.  This should not slow things down very much and does */
/* not really waste any memory. */

/* Implementation functions */

static int RTOp_ROp_num_bounded_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl,  const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  RTOp_value_type        inf_bnd;
  RTOp_index_type        sub_dim;
  const RTOp_value_type  *xl_val;
  ptrdiff_t              xl_val_s;
  const RTOp_value_type  *xu_val;
  ptrdiff_t              xu_val_s;
  RTOp_index_type        num_bounded = 0;
  register RTOp_index_type k;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 2 )
    return RTOp_ERR_INVALID_NUM_VECS;
  assert( vecs );
  if( num_targ_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  if( vecs[0].sub_dim != vecs[1].sub_dim )                 /* Same sizes */
    return RTOp_ERR_INCOMPATIBLE_VECS;

  /* */
  /* Get pointers to data */
  /* */

  /* inf_bnd */
  inf_bnd = *((RTOp_value_type*)obj_data);
  /* sub_dim */
  sub_dim        = vecs[0].sub_dim;
  /* xl */
  xl_val         = vecs[0].values;
  xl_val_s       = vecs[0].values_stride;
  /* xl */
  xu_val         = vecs[1].values;
  xu_val_s       = vecs[1].values_stride;

  /* */
  /* Count the number of bounded variables */
  /* */
  for( k = 0; k < sub_dim; ++k, xl_val += xl_val_s, xu_val += xu_val_s ) {
    if( *xl_val > -inf_bnd || *xu_val < +inf_bnd )
      ++num_bounded;
  }

  /* */
  /* Add this to the result */
  /* */
  *((RTOp_value_type*)targ_obj) += num_bounded;

  return 0; /* success? */
}

/* Virtual function table pointer */
const struct RTOp_RTOp_vtbl_t RTOp_ROp_num_bounded_vtbl =
{
  &RTOp_obj_value_vtbl  /* use simple scalar value type for object instance data */
  ,&RTOp_obj_value_vtbl /* use simple scalar value type for target object */
  ,"ROp_num_bounded"
  ,NULL
  ,RTOp_ROp_num_bounded_apply_op
  ,RTOp_reduct_sum_value
  ,RTOp_get_reduct_sum_value_op
};

/* Class specific functions */

int RTOp_ROp_num_bounded_construct( RTOp_value_type inf_bnd, struct RTOp_RTOp* op )
{
  op->vtbl = &RTOp_ROp_num_bounded_vtbl;
  op->vtbl->obj_data_vtbl->obj_create( NULL, NULL, &op->obj_data );
  *((RTOp_value_type*)op->obj_data) = inf_bnd;
  return 0; /* success? */
}

int RTOp_ROp_num_bounded_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->vtbl      = NULL;
  return 0; /* success? */
}

int RTOp_ROp_num_bounded_set_inf_bnd( RTOp_value_type inf_bnd, struct RTOp_RTOp* op )
{
  *((RTOp_value_type*)op->obj_data) = inf_bnd;
  return 0; /* success? */
}

RTOp_index_type RTOp_ROp_num_bounded_val(RTOp_ReductTarget targ_obj)
{
  return (RTOp_index_type)*((RTOp_value_type*)targ_obj);
}
