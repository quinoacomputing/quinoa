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

// //////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_apply_op_helper.hpp

#ifndef APPLY_OP_HELPER_H
#define APPLY_OP_HELPER_H

#include "AbstractLinAlgPack_Types.hpp"
#include "RTOpPack_RTOpT.hpp"

namespace AbstractLinAlgPack {

/** \brief Validate the inputs to <tt>apply_op()</tt>.
 *
 * Throws an exception if one of the preconditions is not met.
 *
 * ToDo: Finish documentation.
 */
void apply_op_validate_input(
  const char func_name[]
  ,const RTOpPack::RTOp& op
  ,const size_t num_vecs, const Vector* vecs[]
  ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
  );

/** \brief Implements reduction/transformation operators for any serial vectors
 * using just the public vector interface.
 *
 * Note that this function does not validate the input arguments so it is up to
 * the client to do that (i.e. by calling \c apply_op_validate_input()).
 *
 * ToDo: Finish documentation!
 */
void apply_op_serial(
  const RTOpPack::RTOp& op
  ,const size_t num_vecs, const Vector* vecs[]
  ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
  ,RTOpPack::ReductTarget *reduct_obj
  ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
  );

} // end namespace AbstractLinAlgPack

#endif // APPLY_OP_HELPER_H
