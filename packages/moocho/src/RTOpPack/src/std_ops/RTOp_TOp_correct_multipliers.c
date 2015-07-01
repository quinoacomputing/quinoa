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
/* ///////////////////////////////////////////// */
/* RTOp_TOp_Correct_Multipliers.c */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/1/2002 at 18:11 */
/* */

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )

#include "RTOp_TOp_correct_multipliers.h"
#include "RTOp_obj_value_index_vtbl.h"  /* vtbl for operator object instance data */



/* Implementation functions for RTOp_RTOp */

static int RTOp_TOp_Correct_Multipliers_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget reduct_obj )
{
  /* */
  /* Declare local variables */
  /* */

  /* Access to the operator object instance data */
  struct RTOp_value_index_type *data_cntr = (struct RTOp_value_index_type*)obj_data;
  RTOp_value_type *inf_bound_limit = &data_cntr->value;
  RTOp_index_type *lower_or_upper = &data_cntr->index;
  /* Vector data */
  RTOp_index_type           sub_dim;
  /* z0 */
  RTOp_value_type           *z0_val;
  ptrdiff_t                 z0_val_s;
  /* v0 */
  const RTOp_value_type     *v0_val;
  ptrdiff_t                 v0_val_s;

  /* Automatic temporary variables */
  register RTOp_index_type  k;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 1 || ( num_vecs && vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 1 || ( num_targ_vecs && targ_vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  if( /* Validate sub_dim */
    targ_vecs[0].sub_dim != vecs[0].sub_dim
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


  /* */
  /* Apply the operator: */
  /* */
  for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s, z0_val += z0_val_s )
    {
    /* Element-wise transformation */
    if ( (*z0_val) < 0 )
      {
      /* should only get here if the value is slightly negative */
      assert(*z0_val > -1e-50);
      (*z0_val) = 0.0;
      }
      else if ( (*lower_or_upper) == 0 )
        {
      (*z0_val) = ((*v0_val) <= (*inf_bound_limit)) ? 0.0 : (*z0_val);
      }
      else
        {
      (*z0_val) = ((*v0_val) >= (*inf_bound_limit)) ? 0.0 : (*z0_val);
      }
        }

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_TOp_Correct_Multipliers_vtbl =
{
  &RTOp_obj_value_index_vtbl
  ,&RTOp_obj_null_vtbl
  ,"TOp_Correct_Multipliers"
  ,NULL
  ,RTOp_TOp_Correct_Multipliers_apply_op
  ,NULL
  ,NULL
};

/* Class specific functions */

int RTOp_TOp_Correct_Multipliers_construct( RTOp_value_type inf_bound_limit, RTOp_index_type lower_or_upper,  struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert(op);
#endif
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_TOp_Correct_Multipliers_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  return RTOp_TOp_Correct_Multipliers_init(inf_bound_limit,lower_or_upper,op);
}

int RTOp_TOp_Correct_Multipliers_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0;
}

int RTOp_TOp_Correct_Multipliers_init( RTOp_value_type inf_bound_limit, RTOp_index_type lower_or_upper, struct RTOp_RTOp* op )
{
  struct RTOp_value_index_type *ptr_data_cntr = (struct RTOp_value_index_type*)op->obj_data;
  RTOp_value_type *ptr_inf_bound_limit = &ptr_data_cntr->value;
  RTOp_index_type *ptr_lower_or_upper = &ptr_data_cntr->index;
  *ptr_inf_bound_limit = inf_bound_limit;
  *ptr_lower_or_upper = lower_or_upper;
  return 0;
}

