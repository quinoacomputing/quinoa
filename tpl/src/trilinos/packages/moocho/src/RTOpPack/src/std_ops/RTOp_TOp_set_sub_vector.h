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

#ifndef RTOP_TOP_SET_SUB_VECTOR_H
#define RTOP_TOP_SET_SUB_VECTOR_H

#include "RTOp.h"
#include "RTOp_SparseSubVector.h"

#ifdef __cplusplus
extern "C" {
#endif


/** \file RTOp_TOp_set_sub_vector.h Transforamtion operator for setting a sub-vector in a whole vector.
 *
 * <tt>z[0](sub_vec.global_offset+1,sub_vec.global_offset+sub_vec.sub_dim) <- sub_vec</tt>
 *
 * This operator is only defined to allow one vector argument
 * (<tt>num_targ_vecs == 1</tt>) <tt>z[0]</tt>.
 * Using a reduction operator to set a sub-vector will be reasonably
 * efficient for some types of vector subclasses (e.g. dense and sparse serial vectors
 * and parallel vectors with client in every process) but very slow for others
 * (e.g. out-of-core vectors and parallel vectors where operator object must be scattered
 * to various processes).
 * It would be better for vector subclasses to implement this operation directly
 * but if they don't you can use this operator to set the required sub-vector.
 *
 * This operator class works by allocating an internal sub-vector as the operator object
 * state data.
 *
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_set_sub_vector_vtbl;

/* */
/** Constructor.
 *
 * Note that a copy of sub_vec is not made.  Therefore, the client must
 * not disturb <tt>sub_vec</tt> while this operator object is in use!
 */
int RTOp_TOp_set_sub_vector_construct(
  const struct RTOp_SparseSubVector* sub_vec, struct RTOp_RTOp* op );

/* */
/** Reinitialize the range for the sub-vector to extract.
 *
 * Note that a copy of sub_vec is not made.  Therefore, the client must
 * not disturb <tt>sub_vec</tt> while this operator object is in use!
 */
int RTOp_TOp_set_sub_vector_set_sub_vec(
  const struct RTOp_SparseSubVector* sub_vec, struct RTOp_RTOp* op );

/* */
/** Destructor.
 */
int RTOp_TOp_set_sub_vector_destroy( struct RTOp_RTOp* op );

/*@} */

#ifdef __cplusplus
}
#endif

#endif /* RTOP_TOP_SET_SUB_VECTOR_H */
