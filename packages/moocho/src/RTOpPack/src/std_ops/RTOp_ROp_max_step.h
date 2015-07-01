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

#ifndef RTOP_ROP_MAX_STEP_H
#define RTOP_ROP_MAX_STEP_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @name Reduction operator for finding the maximum step for feasibility.
 *
 * <tt>targ_obj <- { max alpha | v[0] + alpha * v[1] >= beta }</tt>
 *
 * This is a specialized reduction operation that is used in many optimization
 * methods to find the maximum step length \c alpha such that the
 * iterates remain feasible.  It is assumed that <tt>v[0] > beta</tt>
 * so that <tt>alpha == 0.0</tt> is a valid return.  If the step is
 * unrestricted, the function \c RTOp_ROp_max_step_val() will return
 * \c RTOp_ROp_max_step_inf.
 *
 * This operator is defined to allow exactly two vector arguments
 * (<tt>num_vecs == 2</tt>) \c v[0], \c v[1], and can only handle
 * dense vectors.
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_step_vtbl;

/* Constructor */
int RTOp_ROp_max_step_construct( RTOp_value_type beta, struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_max_step_destroy( struct RTOp_RTOp* op );

/* Reset beta */
int RTOp_ROp_max_step_set_beta( RTOp_value_type beta, struct RTOp_RTOp* op );

/* Value of max step if step is unrestricted. */
extern RTOp_value_type  RTOp_ROp_max_step_inf;

/* Extract the concrete value from a reduction target object. */
RTOp_value_type  RTOp_ROp_max_step_val(RTOp_ReductTarget reduct_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_ROP_MAX_STEP_H */
