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

#ifndef RTOP_ROP_GET_SUB_VECTOR_H
#define RTOP_ROP_GET_SUB_VECTOR_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif


/** \file RTOp_ROp_get_sub_vector.h Reduction operator for extracting a
 * sub-vector out of the whole vector.
 *
 * <tt>sub_vec <- v[0](l,u)</tt>
 *
 * This operator is only defined to allow one vector argument
 * (<tt>num_vecs == 1</tt>) <tt>v[0]</tt>.
 * Using a reduction operator to extract a sub-vector will be reasonably
 * efficient for some types of vector subclasses (e.g. serial
 * vectors) but very slow for others (e.g. out-of-core and parallel vectors).
 * It would be better for vector subclasses to implement this operation directly
 * but if they don't you can use this operator to extract the required sub-vector.
 *
 * This operator class works by allocating a dense sub-vector and then initializing
 * the values to zero.  During the reduction of a single vector chunk, nonzero values
 * are set.  After the global reduction is performed, the vector is compacted (i.e.
 * the zero elements are removed if present) by the ::RTOp_ROp_get_sub_vector_val<tt>(...)</tt>
 * function.
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_get_sub_vector_vtbl;

/* Constructor */
int RTOp_ROp_get_sub_vector_construct(
  RTOp_index_type l, RTOp_index_type u, struct RTOp_RTOp* op );

/* Reinitialize the range for the sub-vector to extract */
int RTOp_ROp_get_sub_vector_set_range(
  RTOp_index_type l, RTOp_index_type u, struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_get_sub_vector_destroy( struct RTOp_RTOp* op );

/* */
/** Get the sub-vector.
 *
 * This function should only be called after the complete reduction
 * operation has been finished.  However, <tt>RTOp_reduct_obj_free(...)</tt>
 * should still be called to delete the <tt>reduct_obj</tt> once it is
 * finished being used.
 *
 * @param  reduct_obj       [in/out] On input, this object must at least
 *                          be initialized and may or may not have been through
 *                          any reduction operations.
 *
 * @return  The returned sub-vector object does not own any memory.
 * technically, reduct_obj owns the memory.  However, to transfer
 * memory over the client, call <tt>free(reduct_obj)</tt> and then set
 * <tt>reduct_obj == RTOp_REDUCT_OBJ_NULL</tt> to avoid future mistakes.
 */
struct RTOp_SubVector RTOp_ROp_get_sub_vector_val(
  RTOp_ReductTarget  reduct_obj
  );

/*@} */

#ifdef __cplusplus
}
#endif

#endif /* RTOP_ROP_GET_SUB_VECTOR_H */
