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

/* Define this macro if you want to print profiling results with clock() */
#define RTOP_TO_MPI_SHOW_TIMES 1

#ifndef RTOP_TO_MPI_H
#define RTOP_TO_MPI_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \defgroup RTOpToMPI_grp Helper fuctions for applying
 * reduction/transformation operations in C with MPI.
 */
/*@{ */

/* */
/** Initialize MPI compatible type signature arrays for
 * reduction/transformation operator object instance data and
 * reduction target object data.
 *
 * @param num_values    [in] Number of value members
 * @param num_indexes   [in] Number of index members
 * @param num_chars     [in] Number of character members
 * @param num_entries   [out] Number of entries in output arrays set
 * @param block_lengths [out] array (length<tt> >= RTOp_NUM_DATA_TYPES</tt>)
 * @param displacements [out] array (length<tt> >= RTOp_NUM_DATA_TYPES</tt>)
 * @param datatypes     [out] array (length<tt> >= RTOp_NUM_DATA_TYPES</tt>)
 *
 * See the MPI function <tt>MPI_Type_struct(...)</tt> for a discription of these arrays.
 */
void RTOp_MPI_type_signature(
	const int num_values
	,const int num_indexes
	,const int num_chars
	,int* num_entries
	,int block_lengths[]
	,MPI_Aint displacements[]
	,MPI_Datatype datatypes[]
	);

/* */
/** Fill a compacted representation for a reduction object.
 *
 * @param op          [in] RTOp operator object
 * @param reduct_obj  [in] Reduction object
 * @param num_values  [in] Number of value members
 * @param num_indexes [in] Number of index members
 * @param num_chars   [in] Number of character members
 * @param reduct_obj_ext
 *                    [out] (size = <tt>sizeof(RTOp_value_type)*(3 + num_values)
 *                    + sizeof(RTOp_index_type)*num_indexes + sizeof(RTOp_char_type)*num_chars)</tt>.
 *                    On output, this memory will be set according to as described
 *                    in \c RTOp_reduct_op_func_ptr_t.
 *
 * @return Returns 0 if successful otherwise returns an error code.
 */
int RTOp_extract_reduct_obj_ext_state(
	const struct RTOp_RTOp*   op
	,RTOp_ReductTarget        reduct_obj
	,int                      num_values
	,int                      num_indexes
	,int                      num_chars
	,void*                    reduct_obj_ext
	);

/* */
/** Copy from a compacted representation for a reduction object.
 *
 * @param op          [in] RTOp operator object
 * @param reduct_obj_ext
 *                    [in] (size = <tt>sizeof(RTOp_value_type)*(3 + num_values)
 *                    + sizeof(RTOp_index_type)*num_indexes + sizeof(RTOp_char_type)*num_chars)</tt>
 *                    On input, this memory should be set according to as described
 *                    in <tt>#RTOp_reduct_op_func_ptr_t</tt>.
 * @param reduct_obj  [in/out] Reduction object set to the values in \c reduct_obj_ext 
 *                    (must have been previously constructed)
 *
 * @return Returns 0 if successful otherwise returns an error code.
 */
int RTOp_load_reduct_obj_ext_state(
	const struct RTOp_RTOp*   op
	,const void*              reduct_obj_ext
	,RTOp_ReductTarget        reduct_obj
	);

/* */
/** Apply a reduction operation over a set of local sub-vectors using MPI.
 *
 *	@param	comm
 *				[in] MPI communicator
 *	@param	op	[in] reduction/transformation operator object
 *	@param	root_rank
 *				[in] See below
 *	@param	num_cols
 *				[in] The number of columns for each vector
 *	@param	num_vecs
 *				[in] See <tt>RTOpPack::RTOp::apply_op()</tt>
 *	@param	sub_vecs
 *				[in] Array (size <tt>num_targ_vecs * num_cols</tt>)
 *				of nonmutable subvectors.  The vectors for each column kc
 *              begin at sub_vecs+kc*num_vecs where kc=0...num_cols-1
 *              Can be \c NULL if there are no local vector elements.
 *	@param	num_targ_vecs
 *				[in] See <tt>%RTOpPack::RTOp::apply_op()</tt>
 *	@param	sub_targ_vecs
 *				[in] Array (size <tt>num_targ_vecs * num_cols</tt>)
 *				of mutable subvectors.  The vectors for each column kc
 *              begin at sub_vecs+kc*num_targ_vecs where kc=0...num_cols-1
 *              Can be \c NULL if there are no local vector elements.
 *	@param	reduct_objs
 *				[in/out] Array (size <tt>num_cols</tt>) See below.
 *				If <tt>reduct_objs != NULL</tt>
 *				then on output, <tt>reduct_objs[i]</tt> will contain the reduction target
 *				over all the processes along with the reduction on input if
 *              <tt>reduct_objs[i]</tt> has already been through one or more reductions
 *              already and not reinitialized.
 *
 * @return Returns 0 if successful otherwise returns an error code.
 * This function encapsulates a lot of the details of using MPI to perform
 * reduction operations.  For this function to work properly, each
 * processor must contain a single sub-vector for each vector
 * argument and every process (in the communicator) must call this function.
 * Each process's sub-vector for each vector argument
 * may be a different size and even dense on some processes
 * and sparse in others.  It is also allowed for a process to have
 * empty sub-vectors, but still needs to participate in the collective
 * operation and may want the value of the reduction object returned.
 * The main responsibly of the client is to setup the <tt>sub_vecs[]</tt>
 * and <tt>sub_targ_vecs[]</tt> arguments and to deside what type of
 * reduction operation to perform (\c MPI_Allreduce() or MPI_Reduce()).
 *
 * Let <tt>rank</tt> be the value returned from <tt>MPI_Comm_rank(comm,&rank)</tt> in
 * this process.  The expected arguments to this function depend
 * on the argument <tt>root_rank</tt> and <tt>rank</tt> as follows:
 *
 * <ul>
 *	<li> <tt>root_rank < 0</tt> :  In this case we are performing a
 *		<tt>MPI_Allreduce(...)</tt> reduction operation over all of
 *		the processes with the results collected in all of
 *		the processes.  In this case <tt>reduct_obj!=RTOp_REDUCT_OBJ_NULL</tt>
 *		must be true in all of the processes.
 *		<ul>
 *		<li> <tt>root_rank < 0</tt> in all processes
 *		<li> <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt> in all processes
 *           if a reduction is being performed.
 *		</ul>
 *	<li> <tt>root_rank >= 0</tt> :  In this case we are performing a
 *		<tt>MPI_Reduce(...)</tt> reduction operation over all of the processes
 *		with the result only collected in the process with rank
 *		<tt>root_rank</tt>.  Here all of the processes must pass the same
 *		value for <tt>root_rank</tt>.  The reduction target object is only
 *		passed in for the processes with <tt>rank == root_rank</tt>.
 *		<ul>
 *		<li> <tt>root_rank</tt> same in all of the processes
 *		<li> If <tt>rank == root_rank</tt> then <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt>
 *		<li> If <tt>rank != root_rank</tt> then <tt>reduct_obj == RTOp_REDUCT_OBJ_NULL</tt>
 *		</ul>
 *	</ul>
 *
 * Note that if a reduction operation is not performed then no synchronization of the the processes
 * are performed.  This is the fastest behavior and should be fine for a single thread
 * but may cause problems with multiple threads.
 */
int  RTOp_MPI_apply_op(
	MPI_Comm comm, const struct RTOp_RTOp* op, int root_rank
	,const int num_cols
	,const int num_vecs, const struct RTOp_SubVector sub_vecs[]
	,const int num_targ_vecs, const struct RTOp_MutableSubVector sub_targ_vecs[]
	,RTOp_ReductTarget reduct_objs[]
	);

/*@} */

#ifdef RTOP_TO_MPI_SHOW_TIMES
	/* For use in profiling only */
	extern int RTOp_MPI_apply_op_print_timings;
#endif


#ifdef __cplusplus
}
#endif

#endif /* RTOP_TO_MPI_H */
