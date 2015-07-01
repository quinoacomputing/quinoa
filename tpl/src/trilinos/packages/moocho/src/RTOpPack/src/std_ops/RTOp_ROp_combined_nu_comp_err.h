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

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 6/27/2002 at 15:7 */
/* */

#ifndef RTOp_ROp_combined_nu_comp_err_H
#define RTOp_ROp_combined_nu_comp_err_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_combined_nu_comp_err.h
 *
 \verbatim

 element-wise reduction      : comp_err = max(comp_err, v0(i)*(v3(i)-v1(i), -v0(i)*(v1(i)-v2(i))));

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * This operator calculates an estimate of the complementarity error using a combined
 *  nu for the upper and lower bound inequality constraints. At the soln, a positive
 *  nu indicates the upper bound is active, while a negative nu indicates the lower
 *  bound is active. If nu is zero, then neither bound is active.
 *
 *
 \verbatim
      for every i..
       comp_err = max(comp_err, v(i)*(xu(i)-x(i), -v(i)*(x(i)-xl(i))));
\endverbatim
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_combined_nu_comp_err_vtbl;

/* Constructor */
int RTOp_ROp_combined_nu_comp_err_construct(  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_combined_nu_comp_err_destroy( struct RTOp_RTOp* op );


/* Extract the value of the reduction object */
RTOp_value_type RTOp_ROp_combined_nu_comp_err_val(RTOp_ReductTarget reduct_obj);

/*@} */


/**
 *
 \verbatim

 element-wise reduction      : comp_err = max(comp_err, v0(i)*(v1(i)-v2(i)));

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * This operator calculates the comp_err for a single bound vector
 *  (upper or lower) - it is to be used when either one of the bounds
 *  vectors are all infinite (-/+)
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_combined_nu_comp_err_one_only_vtbl;

/* Constructor */
int RTOp_ROp_combined_nu_comp_err_one_only_construct(  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_combined_nu_comp_err_one_only_destroy( struct RTOp_RTOp* op );


/* Extract the value of the reduction object */
RTOp_value_type RTOp_ROp_combined_nu_comp_err_one_only_val(RTOp_ReductTarget reduct_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_ROp_combined_nu_comp_err_H */
