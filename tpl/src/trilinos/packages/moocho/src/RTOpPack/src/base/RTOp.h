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

/* If the macro RTOp_USE_MPI is defined, then these */
/* declarations will be MPI compatible.  If not then */
/* dummy MPI declarations will be used. */
/* */

#ifndef REDUCT_TRANS_VECTOR_OPERATORS_H
#define REDUCT_TRANS_VECTOR_OPERATORS_H

#include <stddef.h>
#include <stdio.h>
#include <assert.h>

#include "RTOp_MPI_config.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef MPIAPI
#define CALL_API MPIAPI
#else
#define CALL_API 
#endif

typedef Teuchos_Ordinal RTOp_index_type;
typedef double RTOp_value_type;
typedef char RTOp_char_type;

/** \defgroup RTOp_grp Interfaces for generalized vector
 * reduction/transformation operators in C.
 *
 * The purpose of these C classes (structs) and functions is
 * to allow users to specify arbitrary reduction/transformation
 * operations on vectors without requiring the vectors to
 * reveal their implementation details.  The design is
 * motivated partly by the "Visitor" patter (Gamma, 1995).
 * Since this set of components is implemented in C it is very
 * type unsafe and the declarations are somewhat verbose
 * to avoid name clashes.  Implementing this in C however makes
 * this set of components accessible to a wider range of users
 * and vector implementations.  However, implementations of the
 * operators in other languages (i.e. C++) is possible.
 * All public declarations associated with
 * this set of components starts with the prefix <tt>RTOp_</tt> which
 * is short for "Reduction/Transformation Operators".
 *
 * This set of interfaces was designed to allow implementation of a
 * distributed parallel application without the explicit knowledge
 * by the application.
 *
 * In the following discussion, <tt>v[k]</tt>, <tt>x</tt>, <tt>y</tt> and <tt>z</tt> are abstract
 * vector objects of dimension <tt>n</tt>.  Users can define operators
 * to perform reduction and/or transformation operations.  Reduction
 * operations must be applied over all of the elements of a vector and
 * therefore requires communication between nodes in a parallel
 * environment but do not change any of the vectors involved.
 * Transformation operations don't require an operation to see an
 * entire vector and therefore usually don't require any communication
 * between nodes in a parallel environment.  The targets of a
 * transformation operation is a set of one or more vectors
 * which are mutated in some way.
 *
 * The idea is that the user may want to perform reduction
 * operations of the form:
 *
 * <tt>op(v[0]...v[*],z[0]...z[*]) -> reduct_obj</tt>
 *
 * where <tt>reduct_obj</tt> is a single object based on a reduction over
 * all the elements of the vector arguments, or transformation
 * operations of the form:
 *
 * <tt>op(v[0](i)...v[*](i),z[0](i)...z[*](i)) -> z[0](i)...z[*](i), for i = 1...n</tt>
 *
 * The tricky part though, is that the <tt>reduct_obj</tt> object of the
 * reduction operation may be more complex
 * than a single scalar value.  For instance, it could be a
 * <tt>double</tt> and an <tt>int</tt> pair such as in the reduction operation:
 *
 * <tt>min{ |x(i)|, i = 1...n } -> [ x(j_min), j_min ]</tt>
 *
 * or it could perform several reductions at once and store
 * several scalar values such as in:
 *
 * <tt>min_max_sum{ x(i), i = 1...n } -> [ x(j_min), j_min, x(j_max), j_max, x_sum ]</tt> 
 *
 * Transformation operations are much simpler to think about
 * and to deal with.  Some off-the-wall examples of transformation
 * operations that this design will support are:
 *
 * <tt>max{ |x(i)|, |y(i)| } + |z(i)| -> z(i), for i = 1...n</tt>
 *
 * <tt>alpha * |z(i)| / x(i) -> z(i), for i = 1...n</tt>
 *
 * <tt>alpha * x(i) * y(i) + beta * z(i) -> z(i), for i = 1...n</tt>
 *
 * Reduction operations present the more difficult technical
 * challenge since they require information gathered from all of the
 * elements to arrive at the final result.  This design allows
 * operator classes to be defined that can simultaneously perform
 * reduction and transformation operations:
  \verbatim
   op(v[0](i)...v[*](i),z[0](i)...z[*](i)) -> z[0](i)...z[*](i),reduct_obj
      , for i = 1...n
  \endverbatim
 * This design is based on a few assumptions about the reduction
 * and transformation operations and vector implementations themselves.
 * First, we will assume that vectors are stored and manipulated as chunks
 * of sub-vectors (of dimension <tt>sub_dim</tt>) where each sub-vector
 * is sufficiently large.  This design supports dense strided
 * sub-vectors (see RTOp_SubVector and RTOp_MutableSubVector)
 * and is relatively flexible.
 *
 * It is strictly the responsibilities of the vector
 * implementations to determine how these operators are applied.
 * For instance, if we are performing a transformation
 * operation of the form:
 *
 * <tt>op( x(i), y(i), z(i) ) -> z(i), for i = 1...n</tt>
 *
 * where <tt>x</tt>, <tt>y</tt>, and <tt>z</tt> are distributed parallel vectors, then we
 * would assume that the elements would be partitioned onto the various
 * processors with the same local elements stored on each processor so
 * as not to require any communication between processors.
 *
 * In order to maintain the simplicity of the design
 * and interfaces, only one real floating-point element type and one
 * index type are supported.
 * The type for the floating-point numbers and the type for the indices
 * should be compatible with Fortran <tt>DOUBLE PRECISION</tt> and <tt>INTEGER</tt> so that
 * these operations can be implemented (the guts anyway) in Fortran or
 * applied to data created in Fortran.  In essence, we want to allow
 * mixed language programming from the start.  To allow
 * for more data types (such as single precision <tt>REAL</tt> or <tt>COMPLEX</tt> )
 * would make these interfaces much more complicated
 * and would be too much of a burden for users and vector implementors
 * to deal with.  This set of datatypes will benefit the largest number
 * of users.
 */
/*@{ */

/** \defgroup RTOp_grp_1 Public declarations, typedefs and misc functions.
 *
 * See \c RTOp_config.h for the definition of platform specific data types
 * \c RTOp_value_type, \c RTOp_index_type and \c RTOp_char_type.
 */
/*@{ */

typedef void*   RTOp_ReductTarget;     /*< The type for reduction target objects. */
#define         RTOp_REDUCT_OBJ_NULL 0 /*< Value given to a a \c NULL reduction target object */
#define         RTOp_NUM_DATA_TYPES 3  /*< Number of primative data types used in \c RTOp */
/* */
/** The prototype for an external reduction operator function (MPI
 * complient).
 *
 * This is a typedef for a function pointer that gets past to
 * <tt>MPI_Reduce(...)</tt> and <tt>MPI_Allreduce(...)</tt>.
 * Therefore, this is the only MPI specific part of this design.
 * Function pointers of this type are returned from \c
 * RTOp_get_reduct_op().  The prototype for this function is:
 *
 * @param array_in    [in] (void *) Array, length <tt>len</tt>
 * @param array_inout [in/out] (void *) Array, length <tt>len</tt>
 * @param len         [in] (int *)
 * @param datatype    [in] (RTOp_Datatype*) Pointer to reduction object datatype (ala MPI)
 * @return void
 *
 * It is very important to note that the reduction objects passed in
 * in <tt>array_in</tt> and <tt>array_inout</tt> are of a special form
 * of the externalized object states.  The arrays of the externalized
 * state passed in and out of \c RTOp_extract_reduct_obj_state() and
 * \c RTOp_load_reduct_obj_state().  The state arrays for a single
 * reduction object must be compacted into a sinlge object \c
 * reduct_obj_ext as follows:
 \verbatim
 RTOp_value_type
   *num_values  = (RTOp_value_type*)reduct_obj_ext,      // Number of elements in values[] 
   *num_indexes = num_values  + sizeof(RTOp_value_type), // Number of elements in indexes[] 
   *num_chars   = num_indexes + sizeof(RTOp_value_type); // Number of elements in chars[] 
 RTOp_value_type
   *values  = num_chars + sizeof(RTOp_value_type);       // Array of num_values values 
 RTOp_index_type
   *indexes = (RTOp_index_type*)(values+num_values);     // Array of num_indexes indexes 
 RTOp_char_type
   *chars   = (RTOp_char_type*)(indexes+num_indexes);    // Array of num_char characters 
 \endverbatim
 * It may seem silly to delcare integer numbers as floating point
 * numbers but the above specification should ensure that the object
 * pointed to by <tt>reduct_obj_ext</tt> will be portable in any
 * heterogeneous environment if we assume that
 * <tt>sizeof(RTOp_value_type) >= sizeof(RTOp_index_type) >=
 * sizeof(RTOp_char_type)</tt> which will be true on most platforms.
 * Arranging the members this way will ensure that none of the
 * individual members will be out of alignment.  It is the
 * responsibility of the vector implementation to create a compacted
 * version of the externalized states of the reduction objects. The
 * only portable way to ensure that the object pointed to by
 * <tt>obj</tt> will be compatible with the above specification is to
 * allocate it as:
 \verbatim
 void *obj = mallac(
               sizeof(RTOp_value_type)*(3 + num_values) +
               sizeof(RTOp_index_type)*num_indexes +
               sizeof(RTOp_char_type)*num_chars
               );
 \endverbatim 
 * Allocating objects in the above way (or by some means equivalent)
 * and explicitly casting the get the individual members may be a
 * little tedious but it is the only way to insure that the object
 * will be layed out properly.  On most platforms with most C
 * compilers however, is may be possible to define structs that will
 * be consistent with the above memory layout of its members but this
 * is not guaranteed by the C standard.
 */

typedef void (CALL_API *RTOp_reduct_op_func_ptr_t) ( void *, void *, int *, RTOp_Datatype * ); 

/** @name Error codes returned from various RTOp functions and interfaces.
 */
/*@{ */
/* */
#define RTOp_ERR_INVALID_USAGE            -1
/* */
#define RTOp_ERR_INVALID_NUM_VECS         -2
/* */
#define RTOp_ERR_INVALID_NUM_TARG_VECS    -3
/* */
#define RTOp_ERR_INCOMPATIBLE_VECS        -5
/* */
#define RTOp_SERVER_INCOMPATIBLE_OPS      -6
/* */
#define RTOp_SERVER_OP_NAME_TOO_LONG      -7
/*@} */
/* */
/** Struct for a non-mutable sub-vector.
 *
 * For a sub-vector <tt>vec</tt>, the corresponding entries
 *	in the global vector <tt>x(j)</tt> (one based) are as follows:
  \verbatim

	x( vec.global_offset + k )
		= vec.values[ vec.value_stride * (k-1) ]

	for k = 1,...,vec.sub_dim
  \endverbatim
 * The stride member <tt>vec.value_stride</tt> may be positive (>0), negative (<0)
 * or even zero (0).  A negative stride <tt>vec.value_stride < 0</tt> allows a
 * reverse traversal of the elements in <tt>vec.values[]</tt>.  A zero stride
 * <tt>vec.value_stride == 0</tt> allows a vector with all the elements the same.
 *
 * To avoid making mistakes in setting the members of this struct use
 * one of the helper functions <tt>RTOp_sub_vector()</tt>,
 * or <tt>RTOp_sub_vector_null()</tt>.
 */
struct RTOp_SubVector {
	/* Offset for the sub-vector into the global vector */
	RTOp_index_type                  global_offset;
	/* Dimension of the sub-vector */
	RTOp_index_type                  sub_dim;
	/* Array (size min{|<tt>value_stride*sub_nz</tt>|,1}) for the values in the vector */
	const RTOp_value_type            *values;
	/* Stride between elements in <tt>values[]</tt> */
	ptrdiff_t                        values_stride;
};
/* */
/** Struct for a mutable sub-vector.
 *
 *	The corresponding entries in the global vector
 * <tt>x(j)</tt> (one based) are as follows:
  \verbatim

	x( vec.global_offset + k )
		= vec.values[ vec.value_stride * (k-1) ]

	for k = 1...vec.sub_dim
  \endverbatim
 * The stride member <tt>vec.value_stride</tt> may be positive (>0), negative (<0)
 * but not zero (0).  A negative stride <tt>vec.value_stride < 0</tt> allows a
 * reverse traversal of the elements in a vector.  A zero stride
 * <tt>vec.value_stride == 0</tt> would allow a vector with all the elements
 * the same and therefore is not a target of a transformation
 * operation.
 *
 * To avoid making mistakes in setting the members of this struct use
 * one of the helper functions <tt>RTOp_mutable_sub_vector()</tt>
 * or <tt>RTOp_mutable_sub_vector_null()</tt>.
 */
struct RTOp_MutableSubVector {
	/* Offset for the sub-vector into the global vector */
	RTOp_index_type                  global_offset;
	/* Dimension of the sub-vector */
	RTOp_index_type                  sub_dim;
	/* Array (size min{|<tt>value_stride*sub_dim</tt>|,1}) for the values in the vector */
	RTOp_value_type	                 *values;
	/* Stride between elements in <tt>values[]</tt> */
	ptrdiff_t                        values_stride;
};
/* */
/** Set the members for a non-mutable sub-vector.
 */
void RTOp_sub_vector(
	RTOp_index_type global_offset, RTOp_index_type sub_dim
	,const RTOp_value_type values[], ptrdiff_t values_stride
	,struct RTOp_SubVector *sub_vec
	);
/* */
/** Initialize a sub-vector argument to null.
 */
void RTOp_sub_vector_null( struct RTOp_SubVector *sub_vec );
/* */
/** Set the members for a mutable dense sub-vector.
 */
void RTOp_mutable_sub_vector(
	RTOp_index_type global_offset, RTOp_index_type sub_dim
	,RTOp_value_type values[], ptrdiff_t values_stride
	,struct RTOp_MutableSubVector *sub_vec
	);
/* */
/** Initialize a sub-vector argument to null.
 */
void RTOp_mutable_sub_vector_null( struct RTOp_MutableSubVector *sub_vec );

/*@} */

/** \defgroup RTOp_grp_2 Reduction/Transformation Operator Interface Functions (virtual).
 *
 * Functions that act virtual with respect to reduction/transformation
 * operators and are used by clients of abstract vectors and by vector
 * implementations to apply these operators.
 *
 * These functions are used as conveniences and call the virtual
 * functions in the <tt>vtbl</tt> of the RTOp_RTOp object and pass in
 * the <tt>obj_data</tt> pointer.  Therefore, these nonmember functions act
 * polymorphically with respect to the operator object.  These
 * functions could be implemented as macros to allow them to
 * be inlined but we all know the problems with macros. The extra
 * function call should not impose too much extra overhead.  Because
 * of the potential for macro inlining, the client should not (and
 * should never need to) take the address of one of these functions.
 *
 * The functions ::RTOp_get_op_type_num_entries<tt>(...)</tt> and
 * ::RTOp_get_reduct_type_num_entries<tt>(...)</tt> are used for externalizing information
 * about the structure of the instance data for reduction/transformation
 * operator objects and for reduction objects.  These functions are needed to
 * externalize the representation of these objects so that these objects can
 * be copied by the client and passed over networks between heterogeneous
 * computers.  This is needed to allow a client/server usage with
 * reduction/transformation operators (see RTOp_Server).
 *
 * To better understand the functions that deal with the opaque reduction target
 * objects, the relationship between reduction/transformation operators and
 * reduction target objects must be clarified.  A reduction object is intimately
 * associated with and is completely owned by an operator object.  In any process,
 * an reduction object can only come into existance by calling the method
 * ::RTOp_reduct_obj_create<tt>(op,&reduct_obj)</tt>.  Once a <tt>reduct_obj</tt> is 
 * created, the memory footprint of the object is set.  If the operator object
 * is later modified in any way (i.e. by ::RTOp_load_op_state<tt>(op,...)</tt>) the
 * the earlier created <tt>reduct_obj</tt> may no longer be compatible with its <tt>op</tt>
 * object.  The implementation of ::RTOp_load_reduct_obj_state<tt>(op,...,reduct_obj)</tt>
 * is not allowed to change the memory footprint of <tt>reduct_obj</tt>.  All of these
 * restrictions are meant to allow for more simplicity in <tt>RTOp</tt> operator
 * implementations.
 *
 */
/*@{ */

struct RTOp_RTOp;

/* */
/** Return the name (as a null-terminated C-style string) of the operator.
 *
 * This name is used to differentate an operator subclass from all
 * other operator subclasses.  This is an important property
 * needed for a client/server and other advanced computing
 * configurations.
 *
 *	@param op        [in] The polymorphic reduction/transformation operator object
 *	@param op_name   [out] Null-terminated string for the name of the operator type
 *
 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_get_op_name(
	const struct RTOp_RTOp* op
	,const char** op_name
	);

/* */
/** Get the number of members of each datatype in the object's externalized state data.
 *
 * See RTOp_obj_type_vtbl_t for a description of this function.
 *
 *	@param op          [in] The polymorphic reduction/transformation operator object
 *	@param num_values  [out] Number of <tt>RTOp_value_type</tt> members
 *	@param num_indexes [out] Number of <tt>RTOp_index_type</tt> members
 *	@param num_chars   [out] Number of <tt>RTOp_char_type</tt> members
 *
 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_get_op_type_num_entries(
	const struct RTOp_RTOp* op
	,int* num_values
	,int* num_indexes
	,int* num_chars
	);

/* */
/** Externalize the state of the operator object to a portable format.
 *
 * This function allows the state of an arbitrary reduction/transformation
 * operator to be transported across a hetergeneous network and have
 * it reconstructed in other processes (using ::RTOp_load_op_state<tt>(...)</tt>).
 *
 * See RTOp_obj_type_vtbl_t for a description of this function.
 * 
 * @return Returns <tt>0</tt> if successful, <tt>!=0</tt> otherwise.
 */
int RTOp_extract_op_state(
	const struct RTOp_RTOp    *op
	,int                      num_values
	,RTOp_value_type          value_data[]
	,int                      num_indexes
	,RTOp_index_type          index_data[]
	,int                      num_chars
	,RTOp_char_type           char_data[]
	);
/* */
/** Load the state of the operator object from a portable format.
 *
 * Note that this function can be called on an uninitilized operator object (i.e.
 * <tt>op->obj_data == NULL</tt>) and in this case, the state data will be dynamicallly
 * allocated in a way that is compatible with the constructors and destructors
 * for the operator class (not given here obviously).  The memory footprint for
 * the operator object may change as a result of this operation even if it has
 * already been initialized.
 *
 * See RTOp_obj_type_vtbl_t for a description of this function.
 *
 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_load_op_state(
	int                       num_values
	,const RTOp_value_type    value_data[]
	,int                      num_indexes
	,const RTOp_index_type    index_data[]
	,int                      num_chars
	,const RTOp_char_type     char_data[]
	,struct RTOp_RTOp*        op
	);
/* */
/** Destroy the state data for this object.
 *
 * @param  op  [in/out] On input, if <tt>op->obj_data != NULL</tt> then this data
 *             will be freed in a way that is compatible with the classes
 *             concrete constructors (not given here of course) and with
 *             ::RTOp_load_op_state<tt>(...)</tt>.  On output, <tt>op->obj_data</tt>
 *             and <tt>op->vtbl</tt> will be set to <tt>NULL</tt>.
 *
 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_free_op( struct RTOp_RTOp* op );
/* */
/** Get the number of members of each datatype in the reduction object.
 *
 * See RTOp_obj_type_vtbl_t for a description of this function.
 *
 *	@param op          [in] The polymorphic reduction/transformation operator object.
 *	@param num_values  [out] Number of <tt>RTOp_value_type</tt> members in the object <tt>op.obj_data</tt>.
 *	@param num_indexes [out] Number of <tt>RTOp_index_type</tt> members in the object <tt>op.obj_data</tt>.
 *	@param num_chars   [out] Number of <tt>RTOp_char_type</tt> members in the object <tt>op.obj_data</tt>.
 *
 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_get_reduct_type_num_entries(
	const struct RTOp_RTOp     *op
	,int                       *num_values
	,int                       *num_indexes
	,int                       *num_chars
	);
/* */
/** Allocate and initialize the reduction object that will be used
 * in the reduction operations.
 *
 * If ::RTOp_get_reduct_type_num_entries<tt>(...)</tt> returns <tt>num_values == 0</tt>
 * , <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then this function should not
 * be called and may be an error if attempted.
 *
 * @param	op	[in] The polymorphic reduction/transformation operator object
 * @param	reduct_obj
 *				[out] On output <tt>*reduct_obj</tt> contains the pointer to
 *				the allocated target object.
 *				Also, <tt>*reduct_obj</tt> will be initialized
 *				ready for use in the reduction operations.
 *				If <tt>*reduct_obj</tt> contains the pointer
 *				to an already allocated object on input, it will not
 *				be freed (see ::RTOp_reduct_obj_free<tt>(...)</tt>).
 *
 * @return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_reduct_obj_create(
	const struct RTOp_RTOp   *op
	,RTOp_ReductTarget       *reduct_obj
	);
/* */
/** Reinitialize an already allocated target object.
 *
 * If ::RTOp_get_reduct_type_num_entries<tt>(...)</tt> returns <tt>num_values == 0</tt>
 * , <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then this function should not
 * be called and may be an error if attempted.
 *
 * @param	op	[in] The reduction/transformation operator object.
 *				This must be the same object that was used in the call to
 *				::RTOp_reduct_obj_create<tt>(op,reduct_obj)</tt>
 * @param	reduct_obj
 *				[out] On output <tt>reduct_obj</tt> will be reinitialized
 *				ready for use in reduction operations.
 *				This object must have been created by the
 *				::RTOp_reduct_obj_create<tt>(op,reduct_obj)</tt>
 *             function first.
 *
 * @return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_reduct_obj_reinit(
	const struct RTOp_RTOp    *op
	,RTOp_ReductTarget        reduct_obj
	);
/* */
/** Free a target object that was previously allocated.
 *
 * If ::RTOp_get_reduct_type_num_entries<tt>(...)</tt> returns <tt>num_values == 0</tt>
 * , <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then this function should not
 * be called and may be an error if attempted.
 *
 * @param	op	[in] The reduction/transformation operator object.
 *				This must be the same object that was used in the call to
 *				::RTOp_reduct_obj_create<tt>(op,reduct_obj)</tt>
 * @param	reduct_obj
 *				[in/out] On input <tt>*reduct_obj</tt> is the pointer to an
 *				allocated target object.  It is allowed that
 *				<tt>*reduct_obj == RTOp_REDUCT_OBJ_NULL</tt>
 *				on input and if so then nothing happens.
 *				This object is then freed and then
 *				on output <tt>*reduct_obj</tt> will be set to <tt>RTOp_REDUCT_OBJ_NULL</tt>
 *
 * @return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_reduct_obj_free( const struct RTOp_RTOp* op
	, RTOp_ReductTarget* reduct_obj );
/* */
/** Externalize the state of the reduction object to a portable format.
 *
 * This allows the state of a reduction object to be transported across
 * a heterogeneous network and also allows the use in MPI global
 * reduction operations.
 *
 * See RTOp_obj_type_vtbl_t for a description of this function.
 *
 * @return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_extract_reduct_obj_state(
	const struct RTOp_RTOp    *op
	,const RTOp_ReductTarget  reduct_obj
	,int                      num_values
	,RTOp_value_type          value_data[]
	,int                      num_indexes
	,RTOp_index_type          index_data[]
	,int                      num_chars
	,RTOp_char_type           char_data[]
	);
/* */
/** Load the state of the reduction object from a portable format.
 *
 * Note that <tt>reduct_obj</tt> must be constructed prior to this and therefore the
 * input data must be compatible with the already constructed <tt>reduct_obj</tt> object.
 *
 * See RTOp_obj_type_vtbl_t for a description of this function.
 *
 * @return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_load_reduct_obj_state(
	const struct RTOp_RTOp   *op
	,int                     num_values
	,const RTOp_value_type   value_data[]
	,int                     num_indexes
	,const RTOp_index_type   index_data[]
	,int                     num_chars
	,const RTOp_char_type    char_data[]
	,RTOp_ReductTarget       reduct_obj
	);
/* */
/** Return if the operator is coordinate invariant.
 *
 * @param  coord_invarient  [out] If <tt>op</tt> is coordinate invarient then
 *                          <tt>*coord_invariant</tt> will be true.
 *
 * @return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_coord_invariant(
	const struct RTOp_RTOp   *op
	,int                     *coord_invariant
	);
/* */
/** <tt>op(sub_vecs[],targ_sub_vecs[]),reduct_obj) -> targ_sub_vecs[],reduct_obj</tt>.
 *
 * This is the bread and butter of the whole design.  Through this method, a
 * vector implementation applies a reduction/transformation operator to a
 * set of sub-vectors.
 *
 * Preconditions:<ul>
 *	<li> <tt>num_vecs > 0 || num_targ_vecs > 0</tt>
 *	<li> <tt>num_vecs > 0 || sub_vecs == NULL</tt>
 *	<li> <tt>num_targ_vecs > 0 || targ_sub_vecs == NULL</tt>
 *	<li> [<tt>num_vecs > 0</tt>] <tt>global_offset == sub_vecs[k].global_offset</tt>
 *        , for <tt>k = 1,...,num_vecs</tt>
 *	<li> [<tt>num_targ_vecs > 0</tt>] <tt>global_offset == targ_sub_vecs[k].global_offset</tt>
 *        , for <tt>k = 1,...,num_targ_vecs</tt>
 *	<li> [<tt>num_vecs > 0</tt>] <tt>sub_dim == sub_vecs[k].sub_dim</tt>
 *       , for <tt>k = 1,...,num_vecs</tt>
 *	<li> [<tt>num_targ_vecs > 0</tt>] <tt>sub_dim == targ_sub_vecs[k].sub_dim</tt>
 *       , for <tt>k = 1,...,num_targ_vecs</tt>
 *	</ul>
 *
 * @param	op	[in] Reduction/transformation operator to apply over the sub-vectors.
 * @param	num_vecs
 *				[in] Number of non-mutable sub-vectors <tt>sub_vec[*]</tt>.
 * @param	sub_vecs
 *				[in] Array (length <tt>num_vecs</tt>) of non-mutable vectors to apply the
 *				operator over.  The ordering of these sub-vectors
 *				<tt>sub_vecs[k], for k = 0...num_vecs-1</tt>, is significant to the <tt>op</tt> object.
 *             If <tt>num_vecs == 0</tt> then <tt>sub_vecs</tt> can be <tt>NULL</tt>.
 * @param	num_targ_vecs
 *				[in] Number of mutable sub-vectors <tt>targ_sub_vec[*]</tt>.
 * @param	targ_sub_vecs
 *				[in] Array (length <tt>num_targ_vecs</tt>) of mutable vectors to apply the
 *				operator over and be mutated.  The ordering of these sub-vectors
 *				<tt>targ_sub_vecs[k], for k = 0...num_targ_vecs-1</tt>, is significant to
 *             the <tt>op</tt> object.  If <tt>num_targ_vecs == 0</tt> then <tt>targ_sub_vecs</tt> can be <tt>NULL</tt>.
 * @param	reduct_obj
 *				[in/out] This reduction object must have been created by
 *				the ::RTOp_reduct_obj_create<tt>(op,reduct_obj)</tt> function and
 *				it may have already passed through one or more other
 *				reduction operations (accumulating the reductions
 *				along the way).  The reduction operation will be:
 *
 *				<tt>op(op(sub_vecs[],targ_sub_vecs[]),reduct_obj) -> reduct_obj</tt>
 *
 *				By allowing an in/out <tt>reduct_obj</tt> and an accumulation
 *				of the reduction, the maximum reuse of memory is achieved.
 *				If <tt>RTOp_reduct_obj_create(op,reduct_obj)</tt> or
 *				::RTOp_reduct_obj_reinit<tt>(op,reduct_obj)</tt> was called
 *				immediately before this function, then <tt>reduct_obj</tt> will
 *				of course only contain the reduction from this operation.
 *             If RTOp_get_reduct_type_num_entries<tt>(...)</tt> returns
 *             <tt>num_values == 0</tt>, <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt>
 *             then <tt>reduct_obj</tt> should be set to <tt>RTOp_REDUCT_OBJ_NULL</tt>
 *             and no reduction will be performed.
 *
 * @return Returns <tt>0</tt> if the operation was successfully executed.
 * If <tt>num_vecs</tt> is incompatible with the underlying operator object then
 * ::RTOp_ERR_INVALID_NUM_VECS is returned and the operation is not performed.
 * If <tt>num_targ_vecs</tt> is incompatible with the underlying operator object then
 * ::RTOp_ERR_INVALID_NUM_TARG_VECS is returned and the operation is not performed.
 * If the sub-vectors are not compatible (i.e. <tt>global_offset</tt> and/or
 * <tt>sub_dim</tt> not the same) then <tt>::RTOp_ERR_INCOMPATIBLE_VECS</tt> is returned.
 */
int RTOp_apply_op(
	const struct RTOp_RTOp                  *op
	,const int                              num_vecs
	,const struct RTOp_SubVector            sub_vecs[]
	,const int                              num_targ_vecs
	,const struct RTOp_MutableSubVector     targ_sub_vecs[]
	,RTOp_ReductTarget                      reduct_obj
	);
/* */
/** <tt>op(in_reduct_obj,inout_reduct_obj) -> inout_reduct_obj</tt>.
 *
 * This function reduces the reduction objects from reduced sub-vectors
 * by the RTOp_apply_op<tt>(op...)</tt> function or those
 * reduced by prior calls to this function.
 *
 * If <tt>reduct_obj  == RTOP_REDUCT_OBJ_NULL</tt> after the return of
 * <tt>RTOp_reduct_obj_create(op,&reduct_obj)</tt>, then this function
 * should not be called and if it is called with arguments that are
 * not <tt>RTOP_REDUCT_OBJ_NULL</tt> then an exception an error value
 * will be returned.
 *
 * @param  op	[in] The reduction/transformation operation used in the calls
 *				to RTOp_apply_op<tt>(op,...)</tt> and prior
 *				calls to this function.
 * @param	in_reduct_obj
 *				[in] A target object from a previous reduction.
 * @param	inout_reduct_obj
 *				[in/out] On input, contains the result from a
 *				previous reduction.  On output, contains the
 *				the reduction of the two target objects.
 *
 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_reduce_reduct_objs(
	const struct RTOp_RTOp     *op
	,RTOp_ReductTarget         in_reduct_obj
	,RTOp_ReductTarget         inout_reduct_obj
	);
/* */
/** Externalize the reduction operation for intermediate target objects.
 *
 * @param	op	[in] The reduction operation used in the calls
 *				to RTOp_apply_op<tt>(op,...)</tt>.
 * @param	reduct_op_func_ptr
 *				[out] On output, <tt>*reduct_op_func_ptr</tt> will
 *				point to an external reduction function
 *				that can be applied to intermediate reduction
 *				target objects.  This function is MPI
 *				compatible and is designed to be used
 *				in MPI reduction operations but may
 *				be used in other contexts.  Any context
 *				specific data needed to perform this reduction
 *				must be contained in the externalized format of
 *				the target objects used with this externalized
 *				reduction function.  It is allowed for an
 *				operator class to return <tt>*reduct_op_func_ptr == NULL</tt>
 *				in which case the client will just have to make due
 *				without this function.
 *
 * @return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
 */
int RTOp_get_reduct_op(
	const struct RTOp_RTOp       *op
	,RTOp_reduct_op_func_ptr_t   *reduct_op_func_ptr
	);

/*@} */

/** \defgroup RTOp_grp_3 Implementation of Reduction/Transformation Operators.
 *
 * These are the structs that must be filled in, in order
 * to create user defined reduction/transformation operators.
 */
/*@{ */

struct RTOp_RTOp_vtbl_t;

/* */
/** Reduction/transformation operation class (struct).
 *
 * Instantiations of this type are used as polymorphic
 * objects for applying reduction/transformation operations
 * on sub-vectors.
 *
 * Strictly speaking, the class of the object is
 * determined by the virtual function table that
 * <tt>vtbl</tt> points to while the specific object
 * instance data is pointed to by <tt>obj_data</tt>.  This
 * design allows complete polymorphic objects
 * in C.
 */
struct RTOp_RTOp {
	/* Pointer to the object data for an instantiation */
	void                           *obj_data;
	/* Pointer to the virtual function table */
	const struct RTOp_RTOp_vtbl_t  *vtbl;
};

/* */
/** Struct for the virtual function table for RTOp_RTOp.
 *
 * This is the table that the user must fill up in order to
 * implement the functions for a reduction operator class.
 *
 * The virtual functions for dealing with the operator instance data and
 * reduction object data are bundled as seperate virtual function tables
 * themselves.  This it to allow for as much reuse as possible since it
 * is expected that the same data structures will be reused an many
 * different situations.
 */
struct RTOp_RTOp_vtbl_t {
	/* Pointer to the virtual function table for the operator object instance data. */
	const struct RTOp_obj_type_vtbl_t  *obj_data_vtbl;
	/* Pointer to the virtual function table for the manipulation of the reduction object. */
	const struct RTOp_obj_type_vtbl_t  *reduct_vtbl;
	/* */
	/** Pointer to a null-terminated string that contains the name of the operator.
	 */
	const char* op_name;
	/* */
	/** Used to overide the initialization or reinitialization of a reduction object
	 * before it is passed through a series of reductions.
	 *
	 * This function pointer should be made <tt>NULL</tt> if the default initialization performed
	 * by <tt>this->reduct_vtbl->obj_create</tt> and <tt>this->reduct_vtbl->obj_reinit</tt> is
	 * sufficient (which will generally be the case).
	 */
	int (*reduct_obj_reinit)(
		const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
		,RTOp_ReductTarget reduct_obj );
/*	/// Called by <tt>RTOp_coord_invariant()</tt> */
/*	int (*coord_invariant) ( */
/*		const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data */
/*		,int *coord_invariant ); */
	/* Called by <tt>RTOp_apply_op()</tt> */
	int (*apply_op)(
		const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
		,const int num_vecs, const struct RTOp_SubVector sub_vecs[]
		,const int num_targ_vecs, const struct RTOp_MutableSubVector targ_sub_vecs[]
		,RTOp_ReductTarget reduct_obj );
	/* Called by <tt>RTOp_reduce_reduct_objs()</tt> */
	int (*reduce_reduct_objs)(
		const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
		,RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj );
	/* Called by <tt>RTOp_get_reduct_op()</tt> */
	int (*get_reduct_op)(
		const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
		,RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr );
};

/* */
/** Vitual function table for manipulating a simple object and externalizing its structure.
 *
 * The functions pointed to in this table create, initialize and destroy a simple object
 * (some struct) and also the structure of the object.  Functions are also declared
 * for allowing the state of the objects to be externalized into a simple portable
 * format.
 * This virtual function pointer table is used for the object instance data for
 * RTOp_RTOp objects and also for the manipulation of reduction objects.  Since
 * it is expected that the same type of structure can be used and reused in several
 * different contexts, they all share this same virtual function table type.
 *
 * The uses of these functions and the meanings of their arguments may depend on the
 * context of where the vtbl is used.
 */
struct RTOp_obj_type_vtbl_t {
	/* */
	/** Returns the number of entries of each type of member that defines the 
	 * externalized object state.
	 *
	 * This function returns the number of entries of values, indexes and characters
	 * (see above) that defines the externalized object state.
	 *
	 * When this vtbl is used to handle the instance data for an
	 * <tt>RTOp_RTOp</tt> class, the function arguments have the following meaning:
	 *	<ul>
	 * <li> <tt>vtbl</tt>        [in] Virtual function table for the operator's data vtbl.
	 *	<li> <tt>instance_data</tt> [in] This is the instance data for the operator
	 *		object <tt>RTOp_RTOp::obj_data</tt> 
	 *	<li> <tt>num_values</tt>  [out] Number of <tt>RTOp_value_type</tt> members in the object.
	 *	<li> <tt>num_indexes</tt> [out] Number of <tt>RTOp_index_type</tt> members in the object.
	 *	<li> <tt>num_chars</tt>   [out] Number of <tt>RTOp_char_type</tt> members in the object.
	 *	</ul>
	 *
	 * When this vtbl is used to handle the reduction object the function arguments
	 * have the following meaning:
	 *	<ul>
	 * <li> <tt>vtbl</tt>        [in] Virtual function table for the reduction objects vtbl.
	 *	<li> <tt>instance_data</tt> [in] This is the instance data for the operator
	 *		object <tt>RTOp_RTOp::obj_data</tt>
	 *	<li> <tt>num_values</tt>  [out] Number of <tt>RTOp_value_type</tt> members in the reduction object.
	 *	<li> <tt>num_indexes</tt> [out] Number of <tt>RTOp_index_type</tt> members in the reduction object.
	 *	<li> <tt>num_chars</tt>   [out] Number of <tt>RTOp_char_type</tt> members in the reduction object.
	 *	</ul>
	 *
	 * If the object is <tt>NULL</tt> and contains no data, the <tt>num_values = 0</tt>,
	 * <tt>num_indexes = 0</tt> and <tt>num_chars = 0</tt> on return.
	 *
	 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
	 */
	int (*get_obj_type_num_entries)(
		const struct RTOp_obj_type_vtbl_t    *vtbl
		,const void*                         instance_data
		,int*                                num_values
		,int*                                num_indexes
		,int*                                num_chars
		);
	/* */
	/** Create (dynamically) an object of this type.
	 *
	 * When this vtbl is used to handle the instance data for an
	 * <tt>RTOp_RTOp</tt> class, the function arguments have the following meaning:
	 *	<ul>
	 * <li> <tt>vtbl</tt>          [in] Virtual function table for the operator's data vtbl.
	 *	<li> <tt>instance_data</tt> [in] This is ignored
	 *	<li> <tt>obj</tt> [out] <tt>*obj</tt> points to the allocated state data to
	 *     <tt>RTOp_RTOp::obj_data</tt>.
	 *	</ul>
	 *
	 * When this vtbl is used to handle the reduction target object for an
	 * <tt>RTOp_RTOp</tt> class, the function arguments have the following meaning:
	 *	<ul>
	 * <li> <tt>vtbl</tt>          [in] Virtual function table for the reduction objects vtbl.
	 *	<li> <tt>instance_data</tt> [in] This is the instance data from
	 *		<tt>RTOp_RTOp::obj_data</tt>
	 *	<li> <tt>obj</tt> [out] The allocated reduction object for the
	 *		<tt>RTOp_RTOp:</tt> operator object.
	 *	</ul>
	 *
	 * If the object is <tt>NULL</tt> and contains no data, then <tt>*obj == NULL</tt>
	 * on return.
	 *
	 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
	 */
	int (*obj_create)(
		const struct RTOp_obj_type_vtbl_t   *vtbl
		,const void                         *instance_data
		,void                               ** obj
		);
	/* */
	/** Reinitialize an object of this type.
	 *
	 * If <tt>get_obj_type_num_entries(ob,...)</tt> returns <tt>num_values = 0</tt>,
	 * <tt>num_indexes = 0</tt> and <tt>num_chars = 0</tt> then this function does nothing
	 * and <tt>obj</tt> must be set to <tt>NULL</tt> if called.
	 *
	 * When this vtbl is used to handle the instance data for a
	 * <tt>RTOp_RTOp</tt> class, this function is optional.
	 * This function could be used in a reinitialization
	 *	function for the <tt>RTOp_RTOp</tt> object
	 *	to reinitialize the object's instance data <tt>obj_data</tt>.
	 *	Or, the pointer to this function could be left as NULL since the
	 *	client of the operator object can never directly call this
	 *	function.
	 *	<ul>
	 * <li> <tt>vtbl</tt>          [in] Virtual function table for the operator's data vtbl.
	 *	<li> <tt>instance_data</tt> [in] Could be data passed into a reinitilaization
	 *		function or more likely it is just ignored.
	 *	<li> <tt>obj</tt> [in/out] The reinitialized object instance data.
	 *	</ul>
	 *
	 * When this vtbl is used to handle the reduction object for an
	 * <tt>RTOp_RTOp</tt> class, the function arguments have the following meaning:
	 *	<ul>
	 * <li> <tt>vtbl</tt>          [in] Virtual function table for the reduction object's vtbl.
	 *	<li> <tt>instance_data</tt> [in] This is the instance data from <tt>RTOp_RTOp::obj_data</tt>
	 *	<li> <tt>obj</tt> [in/out] The target object in reinitialized so some reasonable
	 * value, ready to participate in a reduction operation.
	 *	</ul>
	 *
	 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
	 */
	int (*obj_reinit)(
		const struct RTOp_obj_type_vtbl_t    *vtbl
		,const void                          *instance_data
		,void                                *obj
		);
	/* */
	/** Destroy an object of this type.
	 *
	 * If <tt>get_obj_type_num_entries(ob,...)</tt> returns <tt>num_values = 0</tt>,
	 * <tt>num_indexes = 0</tt> and <tt>num_chars = 0</tt> then this function does nothing
	 * and <tt>obj</tt> must be set to <tt>NULL</tt> if called.
	 *
	 * When this vtbl is used to handle the instance data for an
	 * <tt>RTOp_RTOp</tt> class, this function is optional.
	 * This function could be used in a destructor
	 *	function for the <tt>RTOp_RTOp</tt> class
	 *	to free the object's instance data <tt>obj_data</tt>.
	 *	Or, the pointer to this function could be left as NULL since the
	 *	client of the operator object can never directly call this
	 *	function.
	 *	<ul>
	 * <li> <tt>vtbl</tt>          [in] Virtual function table for the operator's data vtbl.
	 *	<li> <tt>instance_data</tt> [in] This is completely ignored and is not
	 *		needed for anything.
	 *	<li> <tt>obj</tt> [in/out] The instance data from
	 *		<tt>RTOp_RTOp::obj_data</tt> or <tt>RTOp_TransOp::obj_data</tt> to be
	 *		freed and set to <tt>NULL</tt> on output.
	 *	</ul>
	 *
	 * When this vtbl is used to handle the reduction object for an
	 * <tt>RTOp_RTOp</tt> class, the function arguments have the following meaning:
	 *	<ul>
	 * <li> <tt>vtbl</tt>          [in] Virtual function table for the reduction object's vtbl.
	 *	<li> <tt>instance_data</tt> [in] This is the instance data from <tt>RTOp_RTOp::obj_data</tt>
	 *	<li> <tt>obj</tt> [in/out] The target object to be destroyed and set to
	 *		<tt>RTOp_REDUCT_OBJ_NULL</tt> on output.
	 *	</ul>
	 *
	 *	@return Returns <tt>0</tt> if successful and <tt>!=0</tt> otherwise.
	 */
	int (*obj_free)(
		const struct RTOp_obj_type_vtbl_t   *vtbl
		,const void                         *instance_data
		,void                               **obj
		);
	/* */
	/** Externalize the state of the object to a portable format.
	 *
	 * This function must be called on an initalized object.  It is very
	 * important that these data contain all the information needed
	 * to perform the reduction operation using the externalized MPI
	 * compatilble operator function returned from ::RTOp_get_reduct_op<tt>(...)</tt>
	 * when this is used for a reduction object used in a reduction operation.
	 *
	 * @param  vtbl [in] <tt>this</tt> Virtual function table (or NULL).
	 * @param  instance_data
	 *              [in] Points to the instance data for the RTOp operator object (or NULL).
	 * @param  obj  [in[ Reduction object to be externalized.
	 * @param  num_values
	 *              [in] Value returned from <tt>get_obj_type_num_entries(...)</tt>.
	 * @param  value_data
	 *              [out] Array (size <tt>num_values</tt>) of floating point data.
	 * @param  num_indexes
	 *              [in] Value returned from <tt>get_obj_type_num_entries(...)</tt>.
	 * @param  index_data
	 *              [out] Array (size <tt>num_indexes</tt>) of integral data.
	 * @param  num_chars
	 *              [in] Value returned from <tt>get_obj_type_num_entries(...)</tt>.
	 * @param  char_data
	 *              [out] Array (size <tt>num_chars</tt>) of character data.
	 *
	 * @return Returns <tt>0</tt> if successful, <tt>!=0</tt> otherwise.
	 */
	int (*extract_state)(
		const struct RTOp_obj_type_vtbl_t    *vtbl
		,const void                          *instance_data
		,void                                *obj
		,int                                 num_values
		,RTOp_value_type                     value_data[]
		,int                                 num_indexes
		,RTOp_index_type                     index_data[]
		,int                                 num_chars
		,RTOp_char_type                      char_data[]
		);
	/* */
	/** Load the state of the object from a portable format.
	 *
	 * This function can be called on an already constructed object (<tt>*instance_data != NULL</tt>)
	 * or a NULL object (<tt>*instance_data == NULL</tt>).
	 *
	 * @param  vtbl [in] <tt>this</tt> Virtual function table (or NULL).
	 * @param  instance_data
	 *              [in] Points to the instance data for the RTOp operator object (or NULL).
	 * @param  num_values
	 *              [in] Value returned from <tt>get_obj_type_num_entries(...)</tt>.
	 * @param  value_data
	 *              [in] Array (size <tt>num_values</tt>) of floating point data.
	 * @param  num_indexes
	 *              [in] Value returned from <tt>get_obj_type_num_entries(...)</tt>.
	 * @param  index_data
	 *              [in] Array (size <tt>num_indexes</tt>) of integral data.
	 * @param  num_chars
	 *              [in] Value returned from <tt>get_obj_type_num_entries(...)</tt>.
	 * @param  char_data
	 *              [in] Array (size <tt>num_chars</tt>) of character data.
	 * @param  obj
	 *              [in/out] On input if <tt>*obj != NULL</tt> then it will be assumed
	 *              that this is an initialized object and it will be internally
	 *              type casted into an expected type.  If the current object
	 *              pointed to by <tt>*obj</tt> has insufficient size to load
	 *              the given state or <tt>*obj == NULL</tt> on input then
	 *              <tt>free(*obj)</tt> will be called and then
	 *              <tt>*obj = malloc(...)</tt> will be called to allocate the
	 *              proper size object.  Therefore, on output, the pointer <tt>*obj</tt>
	 *              may be changed from its initial value if need be.
	 *
	 * @return Returns <tt>0</tt> if successful, <tt>!=0</tt> otherwise.
	 */
	int (*load_state)(
		const struct RTOp_obj_type_vtbl_t    *vtbl
		,const void                          *instance_data
		,int                                 num_values
		,const RTOp_value_type               value_data[]
		,int                                 num_indexes
		,const RTOp_index_type               index_data[]
		,int                                 num_chars
		,const RTOp_char_type                char_data[]
		,void                                **obj
		);
};

/*@} */

/** \defgroup RTOp_Server RTOp_Server.
 *
 * Singleton object that aids in the transport and reconstruction of
 * reduction/transformation operators in a distributed environment.
 *
 * This is a singleton object that must be initialized in every process in
 * a distributed client/server application where reduction/transformation
 * operators are used.  This object can be used by a abstract vector
 * implementation to transparently take a reduction/transformation
 * operator object created by the client in one process, and then have
 * this operator object transported over a network and be reconstructed in
 * all of the processes where the operator needs to be applied.
 * This object only supplies the basic functions for setting up the
 * <tt>RTOp_Server</tt> object with operator class names and <tt>vtbl</tt> pointers,
 * querying for an operator class
 * name given its <tt>vtbl</tt> pointer and then reconstructing an operator object
 * given its instance data and the name of its operator class.
 *
 * This object is implemented as a singleton object since only one of these
 * objects needs to be instantiated in each process.
 *
 */
/*@{ */

/* */
/** Function for adding an operator class name and virtual function table pointer.
 *
 * The user or client of an abstract vector is responsible for calling these functions
 * to setup the RTOp_Server object with the names and <tt>vtbl</tt> pointers for each of the
 * reduction/transformation operator classes that will be used with the abstract
 * abstract vectors.  This initialization must take place in each of the processes
 * where these operators are used.
 *
 * Several different unrelated clients of abstract vectors can call these functions
 * for the operators they use.  For example, an optimization algorithm and a nonlinear
 * equation solver algorithm may both call these functions to initialize the
 * <tt>RTOp_Server</tt> object for the reduction/transformation operator classes they
 * use.  It is very likely that many of these operator classes may be the same and
 * this interface allows the same operator class to be added more than once with
 * no harm.  Unrelated client packages may however use the same name for different
 * classes by accident in which case these function will flag the problem
 * when they are added.  Therefore, as long as careful names for operator classes
 * are used, then several different client packages can coexist in the same
 * program.
 *
 * Preconditions:
 * <ul>
 *	<li> <tt>strlen(op_class_name) <= RTOp_MAX_REDUCT_TRANS_OP_CLASS_NAME</tt>
 *	<li> <tt>op_class_vtbl != NULL</tt>
 *	</ul>
 *
 *	@param	op_class_name
 *				[in] Name (null terminated string) of the operator subclass.
 *				Should be unique from all other class names.  This name must not be any
 *				longer than RTOp_MAX_REDUCT_TRANS_OP_CLASS_NAME.
 *	@param	op_class_vtbl
 *				[in] Pointer to the virtual function table for the operator
 *				subclass.  This should also be unique from all other classes
 *				but there is no harm if it is not.
 *
 * @return Returns <tt>0</tt> if the operator was successfully added.
 * Returns <tt>>0</tt> if the name <tt>op_class_name</tt> and pointer <tt>op_class_vtbl</tt>
 * have already been added.  Returns <tt><0</tt> if the pair could not be added
 * due to some problem.  If this function returns ::RTOp_SERVER_INCOMPATIBLE_OPS<tt> < 0</tt>
 * then the name <tt>op_class_name</tt> already exists but the pointer <tt>op_class_vtbl</tt>
 * is not the same.  If this function returns ::RTOp_SERVER_OP_NAME_TOO_LONG<tt> < 0</tt>
 * then <tt>op_class_name</tt> was longer than ::RTOp_MAX_REDUCT_TRANS_OP_CLASS_NAME.
 */
int RTOp_Server_add_op_name_vtbl(
	const char                       op_class_name[]
	,const struct RTOp_RTOp_vtbl_t   *op_class_vtbl
	);

/* */
/** Function for lookup up an operator class name given its <tt>vtbl</tt> pointer.
 *
 * This function is used by a vector implementation to lookup the name of a
 * reduction/transformation operator class given a pointer to its <tt>vtbl</tt>.
 * This is needed in order to communicate the class name, along with the object
 * instance data for an operator to other processes so that the operator
 * object can be reconstructed and applied in these processes.
 *
 *	@param	op_class_vtbl
 *				[in] Pointer to the virtual function table for the operator
 *				subclass.  This name must have been used in a previous call to
 *				<tt>RTOp_Server_add_op_name_vtbl(...)</tt>
 *	@param	op_class_name
 *				[out] Name (null terminated string) of the operator subclass.
 *				This array must contain at least <tt>RTOp_MAX_REDUCT_TRANS_OP_CLASS_NAME+1</tt>
 *				entries that can be set.  Note that there may be several
 *				names associated with a single <tt>op_class_vtbl</tt>.  Which name
 *				is returned is implementation dependent.
 *
 * @return Returns <tt>0</tt> if the operator class <tt>vtbl</tt> was successfully found.
 * Returns <tt><0</tt> if <tt>vtbl</tt> could not be found.
 */
int RTOp_Server_lookup_op_name(
	const struct RTOp_RTOp_vtbl_t   *op_class_vtbl
	,char                           op_class_name[]
	);

/* */
/** Function for constructing an operator given its class name and instance data.
 *
 * This function is used by the vector implementation to reconstruct reduction
 * /transformation operators on remote processors where they need to be
 * applied.
 *
 *	@param  op_class_name
 *	                    [in] Name (null terminated string) of the operator subclass.
 *                     This name must have been used in a previous call to
 *                     <tt>RTOp_Server_add_op_name_vtbl(...)</tt>.
 * @param  num_values  [in] Returned from ::RTOp_extract_op_state<tt>(...)</tt>
 * @param  value_data  [in] Returned from <tt>RTOp_extract_op_state(...)</tt>
 * @param  num_indexes [in] Returned from <tt>RTOp_extract_op_state(...)</tt>
 * @param  index_data  [in] Returned from <tt>RTOp_extract_op_state(...)</tt>
 * @param  num_chars   [in] Returned from <tt>RTOp_extract_op_state(...)</tt>
 * @param  char_data   [in] Returned from <tt>RTOp_extract_op_state(...)</tt>
 * @param  op          [out] On input, <tt>op->obj_data == NULL</tt> must be true.
 *                     Constructed operator object.  This object will be
 *                     constructed using ::RTOp_load_op_state<tt>(...)</tt>.  To deallocate
 *                     this object, use ::RTOp_free_op<tt>(...)</tt>.
 *
 * @return Returns <tt>>=0</tt> if the operator object was successfully initialized.
 * Returns <tt><0</tt> if the name <tt>op_class_name</tt> was not recognized.
 */
int RTOp_Server_construct_op(
	const char               op_class_name[]
	,int                     num_values
	,const RTOp_value_type   value_data[]
	,int                     num_indexes
	,const RTOp_index_type   index_data[]
	,int                     num_chars
	,const RTOp_char_type    char_data[]
	,struct RTOp_RTOp        *op
	);

/* */
/** Dump the contents of the name/vtbl pairs for the reduction/transformation
 * operator classes to a C stream (debugging purposes).
 */
void RTOp_Server_dump( FILE* file );

/*@} */

/*@} */

#ifdef __cplusplus
}
#endif

#endif /* REDUCT_TRANS_VECTOR_OPERATORS_H */
