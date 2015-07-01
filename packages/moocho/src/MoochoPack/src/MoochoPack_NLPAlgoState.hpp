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

#ifndef RSQP_STATE_H
#define RSQP_STATE_H

#include <deque>

#include "MoochoPack_Types.hpp"
#include "IterationPack_IterQuantityAccess.hpp"
#include "IterationPack_AlgorithmState.hpp"
#include "IterationPack_cast_iq.hpp"
#include "IterationPack_IterQuantityAccessContiguous.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_Permutation.hpp"
#include "ConstrainedOptPack_DecompositionSystem.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
//#include "DenseLinAlgPack_IVector.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

/** \defgroup rSQPState_IQ_names_grp NLPAlgoState iteration quantities names.
 *
 * ToDo: Finish documentation!
 */
//@{

// Iteration Info
extern const std::string num_basis_name;
// NLP Problem Info 
extern const std::string x_name;
extern const std::string f_name;
extern const std::string Gf_name;
extern const std::string HL_name;
extern const std::string c_name;
extern const std::string h_name;
extern const std::string Gc_name;
// Constraint Gradient Null Space / Range Space Decomposition Info
extern const std::string Y_name;
extern const std::string Z_name;
extern const std::string R_name;
extern const std::string Uy_name;
extern const std::string Uz_name;
// Search Direction Info
extern const std::string py_name;
extern const std::string Ypy_name;
extern const std::string pz_name;
extern const std::string Zpz_name;
extern const std::string d_name;
// Reduced QP Subproblem Info
extern const std::string rGf_name;
extern const std::string rHL_name;
extern const std::string w_name;
extern const std::string zeta_name;
extern const std::string qp_grad_name;
extern const std::string eta_name;
// Global Convergence Info
extern const std::string alpha_name;
extern const std::string merit_func_nlp_name;
extern const std::string mu_name;
extern const std::string phi_name;
// KKT Info
extern const std::string opt_kkt_err_name;
extern const std::string feas_kkt_err_name;
extern const std::string comp_kkt_err_name;
extern const std::string GL_name;
extern const std::string rGL_name;
extern const std::string lambda_name;
extern const std::string nu_name;

//@}

/** \defgroup rSQPState_Macros_grp Macros for adding IterQuantity objects to NLPAlgoState.
 *
 * These macros make it easy to add and remove iteration quantities of any type
 * to the state class.  These macros can even be used by subclasses of NLPAlgoState
 * (only any <tt>IterationPack::AlgorithmState</tt> subclass) to add
 * iteration quantities.  Since scalars and vectors are so pervasive they have
 * there own special macros that take care of actually instantiating the
 * iteration quantities.
 */
//@{

/** \brief Add class declarations for an arbitrary iteration quantity
 */
#define STATE_IQ_DECL(TYPE,NAME)                                          \
  virtual IterQuantityAccess<TYPE>&       NAME();                       \
  virtual const IterQuantityAccess<TYPE>& NAME() const;                 \
private:                                                                  \
  iq_id_encap NAME ## _iq_id_;                                          \
public:

/** \brief Add class declarations for an index (i.e. index_type) iteration quantity
 */
#define STATE_INDEX_IQ_DECL(NAME)                                    \
    STATE_IQ_DECL(index_type,NAME)                                   \

/** \brief Add class declarations for a scalar (i.e. value_type) iteration quantity
 */
#define STATE_SCALAR_IQ_DECL(NAME)                                   \
    STATE_IQ_DECL(value_type,NAME)                                   \

/** \brief Add class declarations for a VectorMutable iteration quantity.
 */
#define STATE_VECTOR_IQ_DECL(NAME)                                   \
    STATE_IQ_DECL(VectorMutable,NAME)                          \

/** \brief Add class definitions for an arbitrary iteration quantity.
 *
 * This implementation can not instantate the underlying iteration quantity
 * so it is the responsibility of some client to do this.  This implementation
 * just initializes the iq_id for the iteration quantity on the fly and
 * then casts the iteration quantity. 
 */
#define STATE_IQ_DEF(CLASS,TYPE,NAME,NAME_STR)                            \
IterQuantityAccess<TYPE>&                                                 \
CLASS::NAME()                                                             \
{                                                                         \
  update_iq_id( NAME_STR, &NAME ## _iq_id_ );                           \
  return IterationPack::cast_iq<TYPE>(                           \
        *this, NAME ## _iq_id_.iq_id, NAME_STR );                         \
}                                                                         \
const IterQuantityAccess<TYPE>&                                           \
CLASS::NAME() const                                                       \
{                                                                         \
  return const_cast<CLASS*>(this)->NAME();                              \
}

/** \brief Add class definitions for a index_type iteration quantity.
 *
 * This implementation will instantiate the IterQuantity object on the fly.
 */
#define STATE_INDEX_IQ_DEF(CLASS,NAME,NAME_STR)                           \
IterQuantityAccess<index_type>&                                           \
CLASS::NAME()                                                             \
{                                                                         \
  update_index_type_iq_id( NAME_STR, &NAME ## _iq_id_ );                \
  return IterationPack::cast_iq<index_type>(                     \
        *this, NAME ## _iq_id_.iq_id, NAME_STR );                         \
}                                                                         \
const IterQuantityAccess<index_type>&                                     \
CLASS::NAME() const                                                       \
{                                                                         \
  return const_cast<CLASS*>(this)->NAME();                              \
}

/** \brief Add class definitions for a value_type iteration quantity.
 *
 * This implementation will instantiate the IterQuantity object on the fly.
 */
#define STATE_SCALAR_IQ_DEF(CLASS,NAME,NAME_STR)                          \
IterQuantityAccess<value_type>&                                           \
CLASS::NAME()                                                             \
{                                                                         \
  update_value_type_iq_id( NAME_STR, &NAME ## _iq_id_ );                \
  return IterationPack::cast_iq<value_type>(                     \
        *this, NAME ## _iq_id_.iq_id, NAME_STR );                         \
}                                                                         \
const IterQuantityAccess<value_type>&                                     \
CLASS::NAME() const                                                       \
{                                                                         \
  return const_cast<CLASS*>(this)->NAME();                              \
}

/** \brief Add class definitions for a VectorMutable iteration quantity.
 *
 * This implementation will instantiate the IterQuantity object on the fly
 * given the VectorSpace (VEC_SPC).  Note that this VEC_SPC can be any
 * code that will returns a smart pointer to a AbstractFactory<VectorMutable>
 * object from within the class body.  It is best if VEC_SPC is some function
 * that is called on *this for the maximum safety and to avoid strage
 * behavior.
 */
#define STATE_VECTOR_IQ_DEF(CLASS,NAME,NAME_STR,VEC_SPC,VEC_RN)           \
IterQuantityAccess<VectorMutable>&                                  \
CLASS::NAME()                                                             \
{                                                                         \
    update_vector_iq_id( NAME_STR, VEC_SPC, VEC_RN, &NAME ## _iq_id_ );   \
  return IterationPack::cast_iq<VectorMutable>(            \
        *this, NAME ## _iq_id_.iq_id, NAME_STR );                         \
}                                                                         \
const IterQuantityAccess<VectorMutable>&                            \
CLASS::NAME() const                                                       \
{                                                                         \
  return const_cast<CLASS*>(this)->NAME();                          \
}

//@}

/** \brief Reduced space SQP state encapsulation interface.
 *
 * This in an interface to a set of data specific to a reduced space SQP algorithms.
 * The iteration quantites are abstracted within <tt>IterQuantityAccess<></tt> objects.
 * A set of boilerplate macros are used to add the necessary declarations and
 * implemetations of these iteration quantity access functions.  As shown by
 * these macros the access methods are declared virtual so that subclasses can
 * override these methods.  Otherwise, much of these could have been declared
 * inline.
 *
 * The implementation defined in this class uses <tt>IterQuantityAccessContiguous<></tt>
 * for iteration quantities of type <tt>index_type</tt>, <tt>value_type</tt> and
 * <tt>VectorMutable</tt> with a default of one storage location.  The default
 * implementation is able to create the <tt>VectorMutable</tt> iteration quantities
 * by using <tt>VectorSpace</tt> objects that the client sets \c this up with.
 *
 * For all other types of iteration quantities (i.e. <tt>MatrixOp</tt> etc.) the
 * client is responsible for setting the iteration quantity object of type
 * <tt>IterQuantityAccess<></tt>.  The client can also change the type of class
 * used for any iteration quantity by simply calling
 * <tt>AlgorithmState::set_iter_quant(...)</tt>.
 *
 * The number of storage locations for any iteration quantity of type
 * <tt>IterQuantityAccessContiguous<></tt> can be changed by fetching
 * the iteration quantity using the access methods defined here and then
 * using <tt>dynamic_cast<></tt> and calling the
 * <tt>IterQuantityAccessContiguous<>::resize(...)</tt> method.
 *
 * Note that the underlying <tt>AlgorithmState</tt> object will not know
 * about the iteration quantity objects with default implementations
 * until the access functions have been called at least once.
 *
 * ToDo: Finish documentation.
 */
class NLPAlgoState
  : public IterationPack::AlgorithmState // doxygen needs full path
{
public:

  /** @name Public Types */
  //@{

  /// Thrown if an iteration quantity is of an invalid type.
  class InvalidType : public std::logic_error
  {public: InvalidType(const std::string& what_arg) : std::logic_error(what_arg) {}};
  
  /** \brief . */
  typedef Teuchos::RCP<const VectorSpace>    vec_space_ptr_t;

  //@}

protected:

  // /////////////////////////////
  // Protected types.

  /** \brief . */
  struct iq_id_encap {
    iq_id_encap() : iq_id(DOES_NOT_EXIST) {}
    iq_id_type iq_id;
  };

public:

  /** @name Constructors/initializers */
  //@{

  // ToDo: Implement all set_space_xx methods to update factories
  // for all vector iteration quantities.

  /// Set the DecompositionSystem object that all share
  STANDARD_COMPOSITION_MEMBERS( DecompositionSystem, decomp_sys );
  /// Set the VectorSpace of x
  STANDARD_CONST_COMPOSITION_MEMBERS( VectorSpace, space_x );
  /// Set the VectorSpace of c
  STANDARD_CONST_COMPOSITION_MEMBERS( VectorSpace, space_c );
  /** \brief Set the VectorSpace of the range space (py).
   * 
   * Calling this method will cause all of the vector iteration
   * quantity objects set in this space to be updated with this
   * vector space (factory) object.
   */
  void set_space_range (const vec_space_ptr_t& space_range );
  vec_space_ptr_t& get_space_range();
  const vec_space_ptr_t& get_space_range() const;
  const VectorSpace& space_range() const;
  /** \brief Set the VectorSpace of the null space (pz).
   * 
   * Calling this method will cause all of the vector iteration
   * quantity objects set in this space to be updated with this
   * vector space (factory) object.
   */
  void set_space_null (const vec_space_ptr_t& space_null );
  vec_space_ptr_t& get_space_null();
  const vec_space_ptr_t& get_space_null() const;
  const VectorSpace& space_null() const;

  /** \brief Construct
   *
   * Initializes num_basis() == 0
   */
  NLPAlgoState(
    const decomp_sys_ptr_t& decomp_sys   = Teuchos::null
    ,const vec_space_ptr_t& space_x      = Teuchos::null
    ,const vec_space_ptr_t& space_c      = Teuchos::null
    ,const vec_space_ptr_t& space_range  = Teuchos::null
    ,const vec_space_ptr_t& space_null   = Teuchos::null
    );

  /** \brief . */
  virtual ~NLPAlgoState() {}

  //@}

  /** @name Iteration Info */
  //@{

    /// num_basis: Counts basis changes durring the algorithm
  STATE_INDEX_IQ_DECL(num_basis)
  
  //@}

  /** @name NLP Problem Info */
  //@{

  /// x:  The current NLP point
  STATE_VECTOR_IQ_DECL(x)
  /// f:  Objective function value
  STATE_SCALAR_IQ_DECL(f)
  /// Gf:  Gradient of the objective function sorted according to current basis selection ( n x 1 )
  STATE_VECTOR_IQ_DECL(Gf)
  /// HL:  Hessian of the Lagrangian ( n x n 
  STATE_IQ_DECL(MatrixSymOp,HL)
  /// c:  DVector of general nonlinear equality constraints ( m x 1 )
  STATE_VECTOR_IQ_DECL(c)
  /// Gc:  Gradient of equality constraints ('c') matrix ( n x m )
  STATE_IQ_DECL(MatrixOp,Gc)

  //@}

  /** @name Constraint Gradient Null Space / Range Space Decomposition Info */
  //@{

  /// Y:  Range space matrix for Gc ([Y  Z] is non-singular) ( n x r )
  STATE_IQ_DECL(MatrixOp,Y)
  /// Z:  Null space matrix for Gc(equ_decomp)' (Gc(equ_decomp)' * Z) ( n x (n-r) )
  STATE_IQ_DECL(MatrixOp,Z)
  /// R:  Represents the nonsingular matrix Gc(equ_decomp)' * Y ( r x r )
  STATE_IQ_DECL(MatrixOpNonsing,R)
  /// Uy:  Represents Gc(equ_undecomp)' * Y ( (m-r) x r )
  STATE_IQ_DECL(MatrixOp,Uy)
  /// Uz:  Represents Gc(equ_undecomp)' * Z ( (m-r) x (m-r) )
  STATE_IQ_DECL(MatrixOp,Uz)

  //@}

  /** @name Search Direction Info */
  //@{

  /// py:  Range space (dependent) QP solution component ( \c space_range, m x 1 )
  STATE_VECTOR_IQ_DECL(py)
  /// Ypy:  Range space (dependent) contribution to search direction (Ypy = Y * py) ( n x 1 )
  STATE_VECTOR_IQ_DECL(Ypy)
  /// pz:  Null space (independent) QP solution component ( \c space_null, (n-m) x 1 )
  STATE_VECTOR_IQ_DECL(pz)
  /// Zpz:  Null space (independent) contribution to the search direction (Zpz = Z * pz) ( n x 1)
  STATE_VECTOR_IQ_DECL(Zpz)
  /// d:  Search direction (d = Zpz + Ypy) ( n x 1 )
  STATE_VECTOR_IQ_DECL(d)

  //@}

  /** @name QP Subproblem Info */
  //@{

  /// rGf:  Reduced gradient of the objective function ( \c space_null, (n-r) x 1 )
  STATE_VECTOR_IQ_DECL(rGf)
  /// rHL:  Reduced Hessian of the Lagrangian function ( <tt>space_null|space_null</tt>, (n-r) x (n-r) )
  STATE_IQ_DECL(MatrixSymOp,rHL)
  /// w:  QP gradient crossterm correction (Z' * HL * Y * py) ( \c space_null, (n-r) x 1 )
  STATE_VECTOR_IQ_DECL(w)
  /// zeta:  QP crossterm dampening parameter [0, 1]
  STATE_SCALAR_IQ_DECL(zeta)
  /// qp_grad:  QP gradient (qp_grad = rGf + zeta * w) ( (n-m) x 1 )
  STATE_VECTOR_IQ_DECL(qp_grad)
  /// eta:  QP relaxation parameter [0, 1]
  STATE_SCALAR_IQ_DECL(eta)

  //@}

  /** @name Global Convergence Info */
  //@{

  /// alpha:  Line seach parameter
  STATE_SCALAR_IQ_DECL(alpha)
  /// merit_func_nlp: Primary merit function for the NLP
  STATE_IQ_DECL(MeritFuncNLP,merit_func_nlp)
  /// mu:  Merit function penalty parameter
  STATE_SCALAR_IQ_DECL(mu)
  /// phi:  Merit function value
  STATE_SCALAR_IQ_DECL(phi)

  //@}

  /** @name KKT Info */
  //@{

  /// Scaled KKT error for optimality ||rGL||
  STATE_SCALAR_IQ_DECL(opt_kkt_err)
  /// Scaled KKT error for feasibility ||c|| and ||hl <= h <= hu||
  STATE_SCALAR_IQ_DECL(feas_kkt_err)
  /// Scaled KKT error for complementarity (bounds)
  STATE_SCALAR_IQ_DECL(comp_kkt_err)
  /// GL:  Gradient of the Lagrangian ( n x 1 )
  STATE_VECTOR_IQ_DECL(GL)
  /// rGL:  Reduced gradient of the Lagrangian ( (n-m) x 1 )
  STATE_VECTOR_IQ_DECL(rGL)
  /// lambda:  Lagrange multipliers for the equality constraints 'c' ( m x 1 )
  STATE_VECTOR_IQ_DECL(lambda)
  /// nu:  Difference between Lagrange multipiers for the upper and lower bounds ( n x 1 )
  STATE_VECTOR_IQ_DECL(nu)

  //@}

  /** @name Decomposition information */
  //@{

  /// Range of dependent variables [1,r]
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Range1D, var_dep );
  /// Range of independent varaibles [r+1,n]
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Range1D, var_indep );

  /// Range of decomposed equality constraints [1,r]
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Range1D, equ_decomp );
  /// Range of undecomposed equality constraints [r+1,m]
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Range1D, equ_undecomp );

  //@}

  /** @name Basis Pivot Info (variable reduction decompositions only) */
  //@{

  /// Current permutation for variables
  STANDARD_COMPOSITION_MEMBERS( Permutation, P_var_current );
  /// Previous permutation for variables
  STANDARD_COMPOSITION_MEMBERS( Permutation, P_var_last );
  /// Current permutation for equality constraints
  STANDARD_COMPOSITION_MEMBERS( Permutation, P_equ_current );
  /// Previous permutation for equality constraints
  STANDARD_COMPOSITION_MEMBERS( Permutation, P_equ_last );

  //@}

protected:

  enum { NUM_VEC_SPACE_TYPES = 5 };
  enum EVecSpaceType {
    VST_SPACE_X       = 0
    ,VST_SPACE_C      = 1
    ,VST_SPACE_RANGE  = 2
    ,VST_SPACE_NULL   = 3
  };

  // /////////////////////////////
  // Protected member functions

  // These implementations are used to avoid code blot and help in debugging
  // (can't debug macros very well).

  /** \brief . */
  void update_iq_id(
    const std::string&                iq_name
    ,iq_id_encap*                     iq_id
    ) const;
  /** \brief . */
  void update_index_type_iq_id(
    const std::string&                iq_name
    ,iq_id_encap*                     iq_id
    );
  /** \brief . */
  void update_value_type_iq_id(
    const std::string&                iq_name
    ,iq_id_encap*                     iq_id
    );
  /** \brief . */
  void update_vector_iq_id(
    const std::string&                iq_name
    ,const VectorSpace::space_ptr_t&  vec_space
    ,EVecSpaceType                    vec_space_type
    ,iq_id_encap*                     iq_id
    );

private:

  // ////////////////////////////
  // Private types

  typedef std::deque<iq_id_type>  iq_vector_list_t;
  
  // ////////////////////////////
  // Private data member

  vec_space_ptr_t    space_range_;
  vec_space_ptr_t    space_null_;

  iq_vector_list_t   vector_iqs_lists_[NUM_VEC_SPACE_TYPES];

  // ////////////////////////////
  // Private member functions.
  
  // Update the vector factories for all of the iteration quantities
  // in the input list.
  void update_vector_factories(
    EVecSpaceType             vec_space_type
    ,const vec_space_ptr_t&   vec_space
    );

  // not defined and not to be called
  NLPAlgoState(const NLPAlgoState&);
  NLPAlgoState& operator=(const NLPAlgoState&);

};	// end class NLPAlgoState

// ////////////////////////////////////
// Inline members

inline
NLPAlgoState::vec_space_ptr_t& NLPAlgoState::get_space_range()
{	return space_range_ ; }

inline
const NLPAlgoState::vec_space_ptr_t& NLPAlgoState::get_space_range() const
{	return space_range_; }

inline
const VectorSpace& NLPAlgoState::space_range() const
{	return *space_range_; }

inline
NLPAlgoState::vec_space_ptr_t& NLPAlgoState::get_space_null()
{	return space_null_ ; }

inline
const NLPAlgoState::vec_space_ptr_t& NLPAlgoState::get_space_null() const
{	return space_null_; }

inline
const VectorSpace& NLPAlgoState::space_null() const
{	return *space_null_; }

}	// end namespace MoochoPack

#endif	// RSQP_STATE_H
