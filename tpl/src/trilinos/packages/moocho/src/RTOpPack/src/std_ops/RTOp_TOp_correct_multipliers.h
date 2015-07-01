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
/* RTOp_TOp_Correct_Multipliers.h */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/1/2002 at 18:11 */
/* */

#ifndef RTOp_TOp_correct_multipliers_H
#define RTOp_TOp_correct_multipliers_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_correct_multipliers.h
 *
 \verbatim


element-wise transformation:
    if (lower_or_upper == 0) // lower bound
        { z0 = (v0 <= inf_bound_limit) ? 0.0 : z0; }
    else // upper bound
        { z0 = (v0 >= inf_bound_limit) ? 0.0 : z0; }

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * This class sets the corresponding multiplier value
 *  to zero if the bound is equal or outside
 *  inf_bound_limit
 *  lower_or_upper : pass 0 for lower bound check
 *                   pass 1 for upper bound check
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_Correct_Multipliers_vtbl;

/* Constructor */
int RTOp_TOp_Correct_Multipliers_construct( RTOp_value_type inf_bound_limit, RTOp_index_type lower_or_upper,  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_TOp_Correct_Multipliers_destroy( struct RTOp_RTOp* op );

/* Initialize the state of the operator object */
int RTOp_TOp_Correct_Multipliers_init( RTOp_value_type inf_bound_limit, RTOp_index_type lower_or_upper, struct RTOp_RTOp* op );



/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_TOp_Correct_Multipliers_H */
