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

#include "RTOpToMPI.h"

#include <stdlib.h>

#ifdef RTOP_TO_MPI_SHOW_TIMES
#include <time.h>
#include <stdio.h>
#endif

void RTOp_MPI_type_signature(
	const int num_values
	,const int num_indexes
	,const int num_chars
	,int* num_entries
	,int block_lengths[]
	,MPI_Aint displacements[]
	,MPI_Datatype datatypes[]
	)
{
	int k = 0, off = 0;
	/* values */
	block_lengths[k] = 3 + num_values; /* must carry size information */
	displacements[k] = 0;
	datatypes[k]     = RTOpMPI_VALUE_TYPE;
	++k;
	off += (3 + num_values) * sizeof(RTOp_value_type);
	/* indexes */
	if( num_indexes ) {
		block_lengths[k] = num_indexes;
		displacements[k] = off;
		datatypes[k]     = RTOpMPI_INDEX_TYPE;
		++k;
		off += num_indexes * sizeof(RTOp_index_type);
	}
	/* chars */
	if( num_chars ) {
		block_lengths[k] = num_chars;
		displacements[k] = off;
		datatypes[k]     = RTOpMPI_CHAR_TYPE;
		++k;
	}
	*num_entries = k;
}

int RTOp_extract_reduct_obj_ext_state(
	const struct RTOp_RTOp*   op
	,RTOp_ReductTarget        reduct_obj
	,int                      num_values
	,int                      num_indexes
	,int                      num_chars
	,void*                    reduct_obj_ext
	)
{
	int err = 0;
	int
		num_values_off  = 0,
		num_indexes_off = num_values_off  + sizeof(RTOp_value_type),
		num_chars_off   = num_indexes_off + sizeof(RTOp_value_type),
		values_off      = num_chars_off   + sizeof(RTOp_value_type),
		indexes_off     = values_off      + num_values  * sizeof(RTOp_value_type),
		chars_off       = indexes_off     + num_indexes * sizeof(RTOp_index_type);
	*(RTOp_value_type*)((char*)reduct_obj_ext + num_values_off)  = num_values;
	*(RTOp_value_type*)((char*)reduct_obj_ext + num_indexes_off) = num_indexes;
	*(RTOp_value_type*)((char*)reduct_obj_ext + num_chars_off)   = num_chars;
	err = RTOp_extract_reduct_obj_state(
		op,reduct_obj
		,num_values,  (RTOp_value_type*)((char*)reduct_obj_ext + values_off)
		,num_indexes, (RTOp_index_type*)((char*)reduct_obj_ext + indexes_off)
		,num_chars,   (RTOp_char_type*)((char*)reduct_obj_ext  + chars_off)
		);
	return err;
}

int RTOp_load_reduct_obj_ext_state(
	const struct RTOp_RTOp*   op
	,const void*              reduct_obj_ext
	,RTOp_ReductTarget        reduct_obj
	)
{
	int err = 0;
	int
		num_values_off  = 0,
		num_values      =                   *(RTOp_value_type*)((char*)reduct_obj_ext + num_values_off),
		num_indexes_off = num_values_off  + sizeof(RTOp_value_type),
		num_indexes     =                   *(RTOp_value_type*)((char*)reduct_obj_ext + num_indexes_off),
		num_chars_off   = num_indexes_off + sizeof(RTOp_value_type),
		num_chars       =                   *(RTOp_value_type*)((char*)reduct_obj_ext + num_chars_off),
		values_off      = num_chars_off + sizeof(RTOp_value_type),
		indexes_off     = values_off  + num_values  * sizeof(RTOp_value_type),
		chars_off       = indexes_off + num_indexes * sizeof(RTOp_index_type);
	err = RTOp_load_reduct_obj_state(
		op
		,num_values,   (RTOp_value_type*)((char*)reduct_obj_ext + values_off)
		,num_indexes, (RTOp_index_type*)((char*)reduct_obj_ext + indexes_off)
		,num_chars,   (RTOp_char_type*)((char*) reduct_obj_ext + chars_off)
		,reduct_obj
		);
	return err;
}

int RTOp_MPI_apply_op(
	MPI_Comm comm, const struct RTOp_RTOp* op, int root_rank
	,const int num_cols
	,const int num_vecs, const struct RTOp_SubVector sub_vecs[]
	,const int num_targ_vecs, const struct RTOp_MutableSubVector sub_targ_vecs[]
	,RTOp_ReductTarget reduct_objs[]
    )
{
	int err = 0;
	int num_reduct_type_values = 0, num_reduct_type_indexes = 0,
		num_reduct_type_chars = 0,  num_reduct_type_entries = 0;
	RTOp_ReductTarget           *i_reduct_objs = NULL;
	int                         target_type_block_lengths[3];
	MPI_Aint                    target_type_displacements[3];
	RTOp_Datatype               target_type_datatypes[3];
	MPI_Datatype                mpi_reduct_ext_type = MPI_DATATYPE_NULL;
	int                         reduct_obj_ext_size = 0;
	char                        *i_reduct_objs_ext = NULL;
	char                        *i_reduct_objs_tmp = NULL;
	RTOp_reduct_op_func_ptr_t   reduct_op_func_ptr  = NULL;
	MPI_Op                      mpi_op = MPI_OP_NULL;
	int                         rank = -1;
	int                         k;
	int                         kc;
#ifdef RTOP_TO_MPI_SHOW_TIMES
    const double secs_per_tick = ((double)1.0) / CLOCKS_PER_SEC;
	clock_t ticks_start_start, ticks_start=0, ticks_end = 0;
#endif

	if( comm == MPI_COMM_NULL ) {
		/* Just do the reduction on this processor and be done with it! */
		if( sub_vecs || sub_targ_vecs ) {
			for( kc = 0; kc < num_cols; ++kc ) {
				err = RTOp_apply_op(
					op, num_vecs, sub_vecs+kc*num_vecs, num_targ_vecs, sub_targ_vecs+kc*num_targ_vecs
					,reduct_objs ? reduct_objs[kc] : RTOp_REDUCT_OBJ_NULL
					);
				if (err) goto ERR_LABEL;
			}
		}
		return 0; /* Success! */
	}

	/* 
	 * Get the description of the datatype of the target object. 
	 * We need this in order to map it to an MPI user defined 
	 * data type so that the reduction operations can be performed 
	 * over all of the of processors. 
	 * 
	 * Get the number of the entries in the array that describes 
	 * target type's datatypes 
	 */
	err = RTOp_get_reduct_type_num_entries(
		op, &num_reduct_type_values, &num_reduct_type_indexes, &num_reduct_type_chars );
	if(err) goto ERR_LABEL;
	num_reduct_type_entries
		= (num_reduct_type_values ? 1 : 0)
		+ (num_reduct_type_indexes ? 1 : 0)
		+ (num_reduct_type_chars ? 1 : 0);

	if( num_reduct_type_entries ) {
#ifdef RTOP_TO_MPI_SHOW_TIMES
		if(RTOp_MPI_apply_op_print_timings) {
			printf("RTOp_MPI_apply_op(...) : timing various MPI calls and other activities\n");
			ticks_start_start = clock();
		}
#endif
		/*
		 * There is a non-null reduction target object so we need to reduce 
		 * it across processors 
		 * 
		 * Allocate the intermediate target object and perform the reduction for the 
		 * vector elements on this processor. 
		 */
        i_reduct_objs = malloc( sizeof(RTOp_ReductTarget) * num_cols );
		for( kc = 0; kc < num_cols; ++kc ) {
			i_reduct_objs[kc] = RTOp_REDUCT_OBJ_NULL;
			RTOp_reduct_obj_create( op, i_reduct_objs + kc );
		}
#ifdef RTOP_TO_MPI_SHOW_TIMES
		if(RTOp_MPI_apply_op_print_timings) {
			printf("calling RTOp_apply_op(...)");
			ticks_start = clock();
		}
#endif
		if( sub_vecs || sub_targ_vecs ) {
			for( kc = 0; kc < num_cols; ++kc ) {
				err = RTOp_apply_op(
					op, num_vecs, sub_vecs+kc*num_vecs, num_targ_vecs, sub_targ_vecs+kc*num_targ_vecs
					,i_reduct_objs[kc]
					);
				if (err) goto ERR_LABEL;
			}
		}
#ifdef RTOP_TO_MPI_SHOW_TIMES
		if(RTOp_MPI_apply_op_print_timings) {
			ticks_end = clock();
			printf(" : cpu time = %g s\n", (ticks_end-ticks_start)*secs_per_tick );
		}
#endif
		/* Get the arrays that describe the reduction object datatype */
		if( num_reduct_type_values == 0 ) ++num_reduct_type_entries; /* Add the block for the sizes */
		RTOp_MPI_type_signature(
			num_reduct_type_values, num_reduct_type_indexes, num_reduct_type_chars
			,&num_reduct_type_entries
			,target_type_block_lengths, target_type_displacements, target_type_datatypes
			);
		/* Translate this target datatype description to an MPI datatype */
#ifdef RTOP_TO_MPI_SHOW_TIMES
		if(RTOp_MPI_apply_op_print_timings) {
			printf("calling MPI_Type_struct(...)");
			ticks_start = clock();
		}
#endif
		err = MPI_Type_struct( num_reduct_type_entries
							   , target_type_block_lengths, target_type_displacements
							   , target_type_datatypes, &mpi_reduct_ext_type );
#ifdef RTOP_TO_MPI_SHOW_TIMES
		if(RTOp_MPI_apply_op_print_timings) {
			ticks_end = clock();
			printf(" : cpu time = %g s\n", (ticks_end-ticks_start)*secs_per_tick );
		}
#endif
		if(err) goto ERR_LABEL;
#ifdef RTOP_TO_MPI_SHOW_TIMES
		if(RTOp_MPI_apply_op_print_timings) {
			printf("calling MPI_Type_commit(...)");
			ticks_start = clock();
		}
#endif
		err = MPI_Type_commit( &mpi_reduct_ext_type );
#ifdef RTOP_TO_MPI_SHOW_TIMES
		if(RTOp_MPI_apply_op_print_timings) {
			ticks_end = clock();
			printf(" : cpu time = %g s\n", (ticks_end-ticks_start)*secs_per_tick );
		}
#endif
		if(err) goto ERR_LABEL;
		/* */
		/* Create an external contiguous representation for the intermediate reduction object */
		/* that MPI can use.  This is very low level but at least the user */
		/* does not have to do it. */
		/* */
		reduct_obj_ext_size =
			(3 + num_reduct_type_values) * sizeof(RTOp_value_type) +
			num_reduct_type_indexes      * sizeof(RTOp_index_type) +
			num_reduct_type_chars        * sizeof(RTOp_char_type);
		i_reduct_objs_ext = malloc( reduct_obj_ext_size * num_cols );
		for( kc = 0; kc < num_cols; ++kc ) {
			RTOp_extract_reduct_obj_ext_state(
				op,i_reduct_objs[kc],num_reduct_type_values,num_reduct_type_indexes,num_reduct_type_chars
				,i_reduct_objs_ext+kc*reduct_obj_ext_size
				);
		}
		/* */
		/* Now perform a global reduction over all of the processors. */
		/* */
		/* Get the user defined reduction operation for MPI to use */
		RTOp_get_reduct_op( op, &reduct_op_func_ptr );
		if( reduct_op_func_ptr != NULL ) {
			/* */
			/* The operator object has returned an MPI compatible reduction function so */
			/* MPI will be of great help in preforming the global reduction :-). */
			/* */
			/* Create the MPI reduction operator object */
			/* */
#ifdef RTOP_TO_MPI_SHOW_TIMES
			if(RTOp_MPI_apply_op_print_timings) {
				printf("calling MPI_Op_create(...)");
				ticks_start = clock();
			}
#endif
			err = MPI_Op_create(
				reduct_op_func_ptr
				,1                   /* The op is communitive? yes == 1! */
				,&mpi_op
				);
#ifdef RTOP_TO_MPI_SHOW_TIMES
			if(RTOp_MPI_apply_op_print_timings) {
				ticks_end = clock();
				printf(" : cpu time = %g s\n", (ticks_end-ticks_start)*secs_per_tick );
			}
#endif
			if(err) goto ERR_LABEL;
			if( root_rank >= 0 ) {
				MPI_Comm_rank( comm, &rank );
				/* Apply the reduction operation over all of the processors */
				/* and reduce to one target object only on the root processor */
				i_reduct_objs_tmp = malloc( reduct_obj_ext_size * num_cols );
				err = MPI_Reduce(
					i_reduct_objs_ext
					,rank == root_rank ? i_reduct_objs_tmp : NULL  /* I don't know if this is MPI legal? */
					,num_cols, mpi_reduct_ext_type
					,mpi_op, root_rank, comm
					);
				if(err) goto ERR_LABEL;
				if( rank == root_rank ) { /* bit-by-bit copy */
					for( k = 0; k < reduct_obj_ext_size * num_cols; ++k ) i_reduct_objs_ext[k] = i_reduct_objs_tmp[k];
				}
			}
			else {
				/* Apply the reduction operation over all of the processors */
				/* and reduce to one target object on each processor */
				i_reduct_objs_tmp = malloc( reduct_obj_ext_size * num_cols );
#ifdef RTOP_TO_MPI_SHOW_TIMES
				if(RTOp_MPI_apply_op_print_timings) {
					printf("calling MPI_Allreduce(...)");
					ticks_start = clock();
				}
#endif
				err = MPI_Allreduce(
					i_reduct_objs_ext, i_reduct_objs_tmp, num_cols
					,mpi_reduct_ext_type, mpi_op, comm
					);
#ifdef RTOP_TO_MPI_SHOW_TIMES
				if(RTOp_MPI_apply_op_print_timings) {
					ticks_end = clock();
					printf(" : cpu time = %g s\n", (ticks_end-ticks_start)*secs_per_tick );
				}
#endif
				if(err) goto ERR_LABEL;
				for( k = 0; k < reduct_obj_ext_size * num_cols; ++k ) i_reduct_objs_ext[k] = i_reduct_objs_tmp[k];
			}
			/* Load the updated state of the reduction target object and reduce. */
			for( kc = 0; kc < num_cols; ++kc ) {
				RTOp_load_reduct_obj_ext_state( op, i_reduct_objs_ext+kc*reduct_obj_ext_size, i_reduct_objs[kc] );
				RTOp_reduce_reduct_objs( op, i_reduct_objs[kc], reduct_objs[kc] );
			}
		}
		else {
			/* */
			/* We must do without the MPI compatible reduction function :-(  We must */
			/* manually perform the reduction operation ourselves.  Note, this will */
			/* not be as efficient as when MPI would do it but it helps take some */
			/* of the burden off of the developers of operator classes. */
			/* */
			/* What we need to do is to gather all of the intermediate reduction */
			/* objects to the root process and then do the reduction pair-wise. */
			/* */
			assert( reduct_op_func_ptr ); /* ToDo: Implement! */
		}
	}
	else {
		/* */
		/* There is a null reduction target object so we don't need to reduce */
		/* it across processors */
		/* */
		if( sub_vecs || sub_targ_vecs ) {
			for( kc = 0; kc < num_cols; ++kc ) {
				err = RTOp_apply_op(
					op, num_vecs, sub_vecs+kc*num_vecs, num_targ_vecs, sub_targ_vecs+kc*num_targ_vecs
					,RTOp_REDUCT_OBJ_NULL
					);
				if (err) goto ERR_LABEL;
			}
		}
	}

ERR_LABEL:

	/* Clean up dynamically allocated memory! */

	if( i_reduct_objs_tmp != NULL )
		free( i_reduct_objs_tmp );
	if( mpi_op != MPI_OP_NULL )
		MPI_Op_free( &mpi_op );
	if( i_reduct_objs_ext != NULL )
		free( i_reduct_objs_ext );
	if( mpi_reduct_ext_type != MPI_DATATYPE_NULL )
		MPI_Type_free( &mpi_reduct_ext_type );
	if( i_reduct_objs != NULL ) {
		for( kc = 0; kc < num_cols; ++kc )
			RTOp_reduct_obj_free( op, &i_reduct_objs[0] );
		free( i_reduct_objs );
	}

	return err;
}

#ifdef RTOP_TO_MPI_SHOW_TIMES
int RTOp_MPI_apply_op_print_timings = 0;
#endif
