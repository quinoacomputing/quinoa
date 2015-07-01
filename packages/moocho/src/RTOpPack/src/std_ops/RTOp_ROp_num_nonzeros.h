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

#ifndef RTOP_ROP_NUM_NONZEROS_H
#define RTOP_ROP_NUM_NONZEROS_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_num_nonzeros.h Reduction operator for counting the number of nonzero elements
 * in a vector.
 *
 * <tt>targ_obj <- size( { i | v[0](i) != 0.0, i=1...n } )</tt>
 *
 * This operator is only defined to allow one vector argument
 * (<tt>num_vecs == 1</tt>) <tt>v[0]</tt> but it can handle a dense or sparse vector.
 * This operator will not count any explicit elments with <tt>v[0](i) == 0.0</tt>
 * in sparse or dense format.  Therefore, this is the exact number of
 * nonzero elements reguardless of storage format.
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_num_nonzeros_vtbl;

/* Constructor */
int RTOp_ROp_num_nonzeros_construct( struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_num_nonzeros_destroy( struct RTOp_RTOp* op );

/* Extract the number of nonzeros */
RTOp_index_type RTOp_ROp_num_nonzeros_val(RTOp_ReductTarget targ_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOP_ROP_NUM_NONZEROS_H */
