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

#ifndef RTOP_TOP_FORCE_IN_BOUNDS_H
#define RTOP_TOP_FORCE_IN_BOUNDS_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_force_in_bounds.h Force the elements in a vector to be within upper and lower bounds.
 *
 * Force <tt>vec[0](i) <= targ_vec[0](i) <= vec[1](i), for i = 1...n</tt>.
 *
 * This operator is only admits dense vectors.  This transformation operation
 * performs the following  (<tt>apply_op(...)</tt>):
 \verbatim
            / xl(i) : if x(i) < xl(i)
 x(i) =     | x(i)  : if xl(i) <= x(i) <= xu(i)
            \ xu(i) : if x(i) > xu(i)
 where:
    x  = targ_vec[0]
  xl = vec[0]
  xu = vec[1]
 \endverbatim
 *
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_force_in_bounds_vtbl;

/* Constructor */
int RTOp_TOp_force_in_bounds_construct( struct RTOp_RTOp* op );

/* Destructor */
int RTOp_TOp_force_in_bounds_destroy( struct RTOp_RTOp* op );

/*@} */


/** This operator is used by the interior point algorithm to push
 *  initial variables sufficiently inside the bounds since the
 *  algorithm assumes that they are ALWAYS within bounds.
 *
 \verbatim

element-wise transformation:
    xl_sb = min(v0 + rel_push*(v1-v0), v0 + abs_push);
    xu_sb = max(v1 - rel_push*(v1-v0), v1 - abs_push);
    if (xl_sb >= xu_sb)
        { z0 = v0 + (v1-v0)/2.0; }
    else if (z0 < xl_sb)
        { z0 = xl_sb; }
    else if (z0 > xu_sb)
        { z0 = xu_sb; }
    // Otherwise, leave it

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 *   xl_sb = min(xl+relative_bound_push*(xu-xl),
 *               xl + absolute_bound_push)
 *   xu_sb = max(xu-relative_bound_push*(xu-xl),
 *               xu - absolute_bound_push)
 *   if (xl_sb > xu_sb) then
 *      x = (xl + (xu-xl)/2
 *   else if (x < xl_sb) then
 *      x = xl_sb
 *   else if (x > xu_sb) then
 *      x = xu_sb
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_force_in_bounds_buffer_vtbl;

/* Constructor */
int RTOp_TOp_force_in_bounds_buffer_construct( RTOp_value_type rel_push, RTOp_value_type abs_push,  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_TOp_force_in_bounds_buffer_destroy( struct RTOp_RTOp* op );

/* Initialize the state of the operator object */
int RTOp_TOp_force_in_bounds_buffer_init( RTOp_value_type rel_push, RTOp_value_type abs_push, struct RTOp_RTOp* op );

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_TOP_FORCE_IN_BOUNDS_H */
