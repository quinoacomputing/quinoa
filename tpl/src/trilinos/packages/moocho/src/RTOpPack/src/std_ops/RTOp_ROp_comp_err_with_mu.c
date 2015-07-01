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

/* ///////////////////////////////////////////// */
/* RTOp_ROp_comp_err_with_mu.c */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/24/2002 at 23:46 */
/* */

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )

#include "RTOp_ROp_comp_err_with_mu.h"
#include "RTOp_obj_value_value_vtbl.h"  /* vtbl for operator object instance data */
#include "RTOp_reduct_max_value.h"


/* Implementation functions for RTOp_RTOp */

static int RTOp_ROp_comp_err_with_mu_apply_op(
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
  RTOp_value_type *mu = &data_cntr->value1;
  RTOp_value_type *inf_bound = &data_cntr->value2;
  /* Access to the reduction object data */
  RTOp_value_type *comp_err = (RTOp_value_type*)reduct_obj;
  /* Vector data */
  RTOp_index_type           sub_dim;
  /* v0 */
  const RTOp_value_type     *v0_val;
  ptrdiff_t                 v0_val_s;
  /* v1 */
  const RTOp_value_type     *v1_val;
  ptrdiff_t                 v1_val_s;
  /* v2 */
  const RTOp_value_type     *v2_val;
  ptrdiff_t                 v2_val_s;
  /* v3 */
  const RTOp_value_type     *v3_val;
  ptrdiff_t                 v3_val_s;
  /* v4 */
  const RTOp_value_type     *v4_val;
  ptrdiff_t                 v4_val_s;

  /* Automatic temporary variables */
  register RTOp_index_type  k;
  /* Temporary element-wise reduction object */
  RTOp_value_type comp_err_ith;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 5 || ( num_vecs && vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 || ( num_targ_vecs && targ_vecs == NULL ) )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  if( /* Validate sub_dim */
    vecs[1].sub_dim != vecs[0].sub_dim
    || vecs[2].sub_dim != vecs[0].sub_dim
    || vecs[3].sub_dim != vecs[0].sub_dim
    || vecs[4].sub_dim != vecs[0].sub_dim
    )
    return RTOp_ERR_INCOMPATIBLE_VECS;
  assert(obj_data);
  assert(reduct_obj);


  /* */
  /* Get pointers to data */
  /* */
  sub_dim       = vecs[0].sub_dim;
  /* v0 */
  v0_val        = vecs[0].values;
  v0_val_s      = vecs[0].values_stride;
  /* v1 */
  v1_val        = vecs[1].values;
  v1_val_s      = vecs[1].values_stride;
  /* v2 */
  v2_val        = vecs[2].values;
  v2_val_s      = vecs[2].values_stride;
  /* v3 */
  v3_val        = vecs[3].values;
  v3_val_s      = vecs[3].values_stride;
  /* v4 */
  v4_val        = vecs[4].values;
  v4_val_s      = vecs[4].values_stride;


  /* */
  /* Apply the operator: */
  /* */
  for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s, v1_val += v1_val_s, v2_val += v2_val_s, v3_val += v3_val_s, v4_val += v4_val_s )
  {
    /* Element-wise reduction */
  comp_err_ith = 0;
  if ((*v1_val) > -(*inf_bound))
    { comp_err_ith = fabs((*v3_val)*((*v0_val)-(*v1_val))-(*mu)); }

  if ((*v2_val) < (*inf_bound))
    { comp_err_ith = max(comp_err_ith, fabs((*v4_val)*((*v2_val)-(*v0_val))-(*mu))); }

    /* Reduction of intermediates */
    (*comp_err) = max( (*comp_err), comp_err_ith );
  }

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_ROp_comp_err_with_mu_vtbl =
{
  &RTOp_obj_value_value_vtbl
  ,&RTOp_obj_value_vtbl
  ,"ROp_comp_err_with_mu"
  ,NULL
  ,RTOp_ROp_comp_err_with_mu_apply_op
  ,RTOp_reduct_max_value
  ,RTOp_get_reduct_max_value_op
};

/* Class specific functions */

int RTOp_ROp_comp_err_with_mu_construct( RTOp_value_type mu, RTOp_value_type inf_bound,  struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert(op);
#endif
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_ROp_comp_err_with_mu_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  return RTOp_ROp_comp_err_with_mu_init(mu,inf_bound,op);
}

int RTOp_ROp_comp_err_with_mu_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0;
}

int RTOp_ROp_comp_err_with_mu_init( RTOp_value_type mu, RTOp_value_type inf_bound, struct RTOp_RTOp* op )
{
  struct RTOp_value_value_type *ptr_data_cntr = (struct RTOp_value_value_type*)op->obj_data;
  RTOp_value_type *ptr_mu = &ptr_data_cntr->value1;
  RTOp_value_type *ptr_inf_bound = &ptr_data_cntr->value2;
  *ptr_mu = mu;
  *ptr_inf_bound = inf_bound;
  return 0;
}
RTOp_value_type RTOp_ROp_comp_err_with_mu_val(RTOp_ReductTarget reduct_obj)
{
  return *((RTOp_value_type*)reduct_obj);
}

