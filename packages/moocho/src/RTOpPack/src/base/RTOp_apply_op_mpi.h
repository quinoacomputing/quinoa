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

#ifndef RTOP_APPLY_OP_MPI_H
#define RTOP_APPLY_OP_MPI_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/* */
/** Function that implements the guts an <tt>apply_op()</tt> method for dense MPI vectors.
 *
 * @param  comm        [in] MPI communicator
 * @param  global_dim  [in] Number of global elements represented by the MPI vector object that
 *                     called this function.
 * @param  local_dim   [in] Number of vector elements stored on this local processor
 * @param  num_cols    [in] The number of columns in each mulit-vector argument.
 * @param  num_vecs    [in] Number of non-mutable vectors involved in this reduction/transformation
 *                     operation (see below).
 * @param  local_vec_ptrs
 *                     [in] Array (size \c num_vecs) of pointers to the local elements of the
 *                     non-mutable vectors (see below).
 * @param  local_vec_strides
 *                     [in] Array (size \c num_vecs) of strides between vector elements in
 *                     \c local_vec_ptrs[p] (see below).
 * @param  local_vec_leading_dim
 *                     [in] Array (size \c num_vecs) of leading dimmensions between colummns 
 *                     for each vector in \c local_vec_ptrs[p] (see below).  Can be <tt>NULL</tt>
 *                     if <tt>num_cols==0</tt>.
 * @param  num_targ_vecs
 *                    [in] Number of mutable vectors involved in this reduction/transformation
 *                     operation (see below).
 * @param  local_targ_vec_ptrs
 *                     [in/out] Array (size \c num_targ_vecs) of pointers to the local elements of the
 *                     mutable vectors to be transformed (see below).
 * @param  local_targ_vec_strides
 *                     [in] Array (size \c num_targ_vecs) of leading dimmensions between colummns 
 *                     for each vector in \c local_targ_vec_ptrs[p] (see below).
 * @param  local_targ_vec_leading_dim
 *                     [in] Array (size \c num_vecs) of strides between vector elements in
 *                     \c local_vec_ptrs[p] (see below).  Can be <tt>NULL</tt>
 *                     if <tt>num_cols==0</tt>.
 * @param  first_ele   [in] Identifies the first global element in the input parallel vector that
 *                     defines the logical sub-vector that the RTOp operator will be applied to.
 * @param  sub_dim     [in] Identifies the number of elements in the input parallel vector that
 *                     defines the logical sub-vector that the RTOp operator will be applied to.
 *                     If <tt>sub_dim == 0</tt> then all of the remaining global elements will
 *                     be included in the logical vector.
 * @param  global_offset
 *                     [in] Identifies where the sub-vector selected by \c first_ele and \c sub_dim
 *                     exists in the logical sub-vector that the RTOp operator will be applied to.
 * @param  op          [in] Reduction/transformation operator to apply over each sub-vector
 *                     and use to add to the reduction target object <tt>reduct_obj</tt> (if
 *                     <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt>).
 * @param  reduct_objs
 *                     [in/out] Array (sizee num_cols) of target objects of the reduction operation.
 *                     These objects must have been created by the <tt>RTOp_reduct_obj_create(op&reduct_objs[kc])</tt>
 *                     function first, for kc=0...num_cols-1.  The reduction operation will be added to <tt>(reduct_objs[])</tt> if
 *                     <tt>(reduct_obj[kc])</tt> has already been through a reduction.  By allowing the info in
 *                     <tt>(reduct_objs[kc])</tt> to be added to the reduction over all of these vectors, the reduction
 *                     operation can be accumulated over a set of abstract vectors which can be useful for implementing
 *                     composite vectors for instance.  If <tt>RTOp_get_reduct_type_num_entries(op,...)</tt> returns
 *                     <tt>num_values == 0</tt>, <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then
 *                     <tt>reduct_obj</tt> should be set to <tt>NULL</tt> and no reduction will be performed.
 *
 * This function takes care of all of the ugly details that goes on under the hood of using \c RTOp operators
 * in a distrubuted parallel application using MPI.
 *
 * This first set of arguments specifies the MPI environment and defines the parallel vector
 * arguments and their local data.
 *
 * The set of arguments \c first_ele, \c sub_dim and \c global_offset defines the logical sub-vectors
 * that the operator will be applied to.  See the function \c RTOp_parallel_calc_overlap() for a
 * description of what these arguments mean.
 *
 * This last set of arguments passes in the \c RTOp operator object \c op and the reduction target
 * object \c reduct_obj.
 * 
 * ToDo: Finish documentation! 
 */
int RTOp_apply_op_mpi(
	MPI_Comm comm
	,RTOp_index_type global_dim, RTOp_index_type local_sub_dim, RTOp_index_type local_offset
	,const int num_cols
	,const int num_vecs,       const RTOp_value_type*   local_vec_ptrs[],       const ptrdiff_t local_vec_strides[],      const ptrdiff_t local_vec_leading_dim[]
	,const int num_targ_vecs,  RTOp_value_type*         local_targ_vec_ptrs[],  const ptrdiff_t local_targ_vec_strides[], const ptrdiff_t local_targ_vec_leading_dim[]
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	,const struct RTOp_RTOp* op
	,RTOp_ReductTarget reduct_objs[]
	);

#ifdef __cplusplus
}
#endif

#endif /* RTOP_APPLY_OP_MPI_H */
