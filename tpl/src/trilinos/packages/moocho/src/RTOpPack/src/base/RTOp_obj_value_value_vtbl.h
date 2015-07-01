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

#ifndef RTOP_OBJ_VALUE_VALUE_VTBL_H
#define RTOP_OBJ_VALUE_VALUE_VTBL_H

#include "RTOp.h"

typedef double crap;

#ifdef __cplusplus
extern "C" {
#endif

/** Virtual function table for a simple object of a <tt>{RTOp_value_type,RTOp_value_type}</tt> pair.
  *
  * The functions do the following:
  * <ul>
  *	<li> <tt>get_obj_type_num_entries(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>num_values</tt>    [out] Returns 2
  *		<li> <tt>num_indexes</tt>   [out] Returns 0
  *		<li> <tt>num_chars</tt>     [out] Returns 0
  *		</ul>
  *	<li> <tt>obj_create(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>obj</tt>           [out] Points an allocated to a <tt>RTOp_value_value_type</tt>
  *                                 object initialized to <tt>{0.0,0.0}</tt>
  *		</ul>
  *	<li> <tt>obj_reinit(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>obj</tt>           [in/out] <tt>RTOp_value_value_type</tt> object reinitialized to <tt>{0.0,0.0}</tt>
  *		</ul>
  *	<li> <tt>obj_free(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>obj</tt>           [in/out] allocated object is freed and <tt>obj</tt>
  *                                 is set to NULL.
  *		</ul>
  *	<li> <tt>extract_state(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *     <li> <tt>obj</tt>           [in] Allocated <tt>RTOp_value_value_type</tt> object.
  *		<li> <tt>num_values</tt>    [in] Must be 2
  *     <li> <tt>value_data</tt>    [out] <ul><li><tt>value_data[0] = ((RTOp_value_value_type*)obj->value1)</tt></li>
  *                                       <li><tt>value_data[1] = ((RTOp_value_value_type*)obj->value2)</tt></li></ul>
  *		<li> <tt>num_indexes</tt>   [in] Must be 0
  *     <li> <tt>index_data</tt>    [out] Must be NULL
  *		<li> <tt>num_chars</tt>     [in] Must be 0
  *     <li> <tt>char_data</tt >    [out] Must be NULL
  *		</ul>
  *	<li> <tt>load_state(...)</tt>
  *		<ul>
  *		<li> <tt>vtbl</tt>          [in] Ignored
  *		<li> <tt>instance_data</tt> [in] Ignored
  *		<li> <tt>num_values</tt>    [in] Must be 2
  *     <li> <tt>value_data</tt>    [in] Array (size 2)
  *		<li> <tt>num_indexes</tt>   [in] Must be 0
  *     <li> <tt>index_data</tt>    [in] Must be NULL
  *		<li> <tt>num_chars</tt>     [in] Must be 0
  *     <li> <tt>char_data</tt >    [in] Must be NULL
  *     <li> <tt>obj</tt>           [in/out] If <tt>*obj == NULL</tt> then
  *                                 <tt>*obj = malloc(sizeof(RTOp_value_value_type))</tt> and
  *                                 <tt>*(RTOp_value_value_type)*obj = {value_data[0],value_data[1]}</tt>
  *		</ul>
  * </ul>
  */
/*@{ */

/* Virtual function table */
extern const struct RTOp_obj_type_vtbl_t   RTOp_obj_value_value_vtbl;

/* Object type structure for a value, index pair. */
struct RTOp_value_value_type {
	/* */
	RTOp_value_type  value1;
	/* */
	RTOp_value_type  value2;
};

/*@} */

#ifdef __cplusplus
}
#endif

#endif /* RTOP_OBJ_VALUE_VALUE_VTBL_H */
