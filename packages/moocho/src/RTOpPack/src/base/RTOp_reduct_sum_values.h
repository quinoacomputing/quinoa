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

#ifndef RTOP_REDUCT_SUM_VALUES_H
#define RTOP_REDUCT_SUM_VALUES_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @name Definitions of reduction functions for summing a list of scalar values.
 *
 * These functions perform a simple sum of a list of scalar objects
 * as defined by the virtual function table \Ref{RTOp_obj_values_vtbl}.
 */
/*@{ */

/* */
/** Use this function for <tt>reduce_reduct_objs</tt> in the RTOp_RTOp_vtbl_t virtual
 * function table.
 */
int RTOp_reduct_sum_values(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, RTOp_ReductTarget in_targ_obj, RTOp_ReductTarget inout_targ_obj );

/* */
/** Use this function for <tt>get_reduct_op</tt> in the RTOp_RTOp_vtbl_t virtual
 * function table.
 */
int RTOp_get_reduct_sum_values_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr );

/*@} */

#ifdef __cplusplus
}
#endif

#endif /* RTOP_REDUCT_SUM_VALUES_H */
