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

#ifndef ALAP_VECTOR_MUTABLE_HPP
#define ALAP_VECTOR_MUTABLE_HPP

#include "AbstractLinAlgPack_Vector.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface for mutable coordinate vectors {abstract}.
 *
 * Objects of this type can act as a target vector of a transformation operation.
 * Similarly to <tt>Vector</tt> this interface contains very few (only one extra) pure
 * virtual methods that must be overridden.  However, more efficient and more general
 * implementations will choose to override more methods.
 *
 * In addition to being able to create non-mutable (\c const) abstract sub-views of a vector
 * object thorugh the \c Vector interface, this interface allows the creation of
 * mutable (non-<tt>const</tt>) sub-views using \c sub_view().  Also, in addition to being
 * able to extract an explicit non-mutable view of some (small?) sub-set of elements, this
 * interface allows a client to either extract a explicit mutable sub-views using
 * \c get_sub_vector() or to set sub-vectors using \c set_sub_vector(). As much
 * as possible, abstract views should be preferred (i.e. \c sub_view()) over explict views (i.e.
 * get_sub_vector() and set_sub_vector()).
 *
 * There are only two pure virtual methods that a concreate <tt>VectorMutable</tt>
 * subclass must override.  The <tt>space()</tt> and <tt>apply_op()</tt> methods from the <tt>Vector</tt>
 * base class inteface must be defined.
 *
 * The non-mutable (<tt>const</tt>) <tt>sub_view(...)</tt> method from the <tt>Vector</tt>
 * interface has a default implementation defined here that will be adequate for most subclasses.
 */
class VectorMutable : virtual public Vector
{
public:

  /** \brief . */
  using Vector::get_sub_vector;
  /** \brief . */
  using Vector::free_sub_vector;

  /** @name Virtual methods with default implementations */
  //@{

  /** \brief Assign the elements of <tt>this</tt> vector to a scalar.
   *
   * The default implementation of this function uses a transforamtion operator class
   * (see RTOp_TOp_assign_scalar.h) and calls <tt>this->apply_op()</tt>.
   */
  virtual VectorMutable& operator=(value_type alpha);

  /** \brief Assign the elements of a vector to <tt>this</tt>.
   *
   * The default implementation of this function uses a transforamtion operator class
   * (see RTOp_TOp_assign_vectors.h) and calls <tt>this->apply_op()</tt>.
   */
  virtual VectorMutable& operator=(const Vector& v);

  /** \brief Default implementation calls <tt>operator=((const &Vector)v)</tt>.
   */
  virtual VectorMutable& operator=(const VectorMutable& v);

  /** \brief Set a specific element of a vector.
   *
   * Preconditions:<ul>
   * <li> <tt>1 <= i <= this->dim()</tt> (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->get(i) == val</tt>
   * </ul>
   *
   * The default implementation uses a transforamtion operator
   * class (see RTOp_TOp_set_ele.h) and calls <tt>this->apply_op()</tt>.
   *
   * @param  i    [in] Index of the element value to set.
   * @param  val  [in] Value of the element to set.
   */
  virtual void set_ele( index_type i, value_type val );

  /** \brief Create a mutable abstract view of a vector object.
   *
   * This is only a transient view of a sub-vector that is to be immediately used
   * and then released by <tt>RCP<></tt>.  This function is declared as
   * non-constant because the object returned has the capacity to alter <tt>this</tt>
   * object.
   *
   * The compatibility of sub-views goes along with the compatibility of sub-spaces
   * (see VectorSpace).  For example, given the vector objects where
   * <tt>x.space().is_compatible(y.space()) == true</tt> then if
   * <tt>x.space().sub_space(rng1)->is_compatible(*y.space().sub_space(rng2)) == true</tt>
   * then the sub-vector views <tt>*x.sub_view(rng1)</tt> and <tt>*y.sub_view(rng2)</tt>
   * should be compatible and can be combined in vector operations.
   *
   * Preconditions:<ul>
   * <li> <tt>rng.in_range(this->dim()) == true</tt> (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * @param  rng  [in] The range of the elements to extract the sub-vector view.
   * 
   * @return  Returns a smart reference counted pointer to a view of the requested
   * vector elements.  It is allowed for the vector implementation to refuse to
   * create arbitrary views in which case this function will return
   * <tt>return.get() == NULL</tt>. In most applications, only specific views are
   * every required.  The default implementation uses the subclass VectorSubView
   * to represent any arbitrary sub-view but this can be inefficient if the sub-view is very
   * small compared this this full vector space but not necessarily.  Note that the underlying
   * vector <tt>this</tt> is not guarrenteed to show the changes made the sub-view
   * <tt>*return .get()</tt> until the smart reference counted pointer <tt>return</tt> is
   * destroyed.
   */
  virtual vec_mut_ptr_t sub_view( const Range1D& rng );

  /** \brief Inline member function that simply calls <tt>this->sub_view(Range1D(l,u))</tt>.
   */
  vec_mut_ptr_t sub_view( const index_type& l, const index_type& u );

  /** \brief Zeros this vector.
   *
   * Calls <tt>operator=(0.0)</tt>.
   */
  virtual void zero();

  /** \brief Adds a linear combination of another vector to this vector object.
   *
   * Calls <tt>this->apply_op()</tt> with an operator class
   * (see RTOp_TOp_axpy.h).
   */
  virtual void axpy( value_type alpha, const Vector& x );

  /** \brief Get a mutable explicit view of a sub-vector.
   *
   * This is only a transient view of a sub-vector that is to be immediately used
   * and then released with a call to \c release_sub_vector().
   *
   * Note that calling this operation might require some internal
   * allocations and temporary memory.  Therefore, it is critical
   * that <tt>this->release_sub_vector(sub_vec)</tt> is called to
   * clean up memory and avoid memory leaks after the sub-vector
   * is used.
   *
   * If <tt>this->get_sub_vector(...,sub_vec)</tt> was previously
   * called on <tt>sub_vec</tt> then it may be possible to reuse this
   * memory if it is sufficiently sized.  The user is
   * encouraged to make multiple calls to <tt>this->get_sub_vector(...,sub_vec)</tt>
   * before <tt>this->release_sub_vector(sub_vec)</tt> to finally
   * clean up all of the memory.  Of course the same <tt>sub_vec</tt> object must be
   * passed to the same vector object for this to work correctly.
   *
   * Changes to the underlying sub-vector are not guarrenteed to become permanent
   * until <tt>this->get_sub_vector(...,sub_vec)</tt> or <tt>this->commit_sub_vector(sub_vec)</tt>
   * is called.
   *
   * Preconditions:<ul>
   * <li> [<tt>!rng.full_range()</tt>] <tt>(rng.ubound() <= this->dim()) == true</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * This method has a default implementation based on a vector reduction operator
   * class (see RTOp_ROp_get_sub_vector.h) and calls <tt>apply_op()</tt>.
   * Note that the footprint of the reduction object (both internal and external state)
   * will be O(<tt>rng.size()</tt>).  For serial applications this is faily adequate and will
   * not be a major performance penalty.  For parallel applications, this will be
   * a terrible implementation and must be overridden if <tt>rng.size()</tt> is large at all.
   * If a subclass does override this method, it must also override <tt>release_sub_vector()</tt>
   * which has a default implementation which is a companion to this method's default
   * implementation.
   *
   * @param  rng      [in] The range of the elements to extract the sub-vector view.
   * @param  sub_vec  [in/out] Mutable view of the sub-vector.  Prior to the
   *                  first call <tt>RTOp_mutable_sub_vector_null(sub_vec)</tt> must
   *                  have been called for the correct behavior.  Technically
   *                  <tt>*sub_vec</tt> owns the memory but this memory can be freed
   *                  only by calling <tt>this->commit_sub_vector(sub_vec)</tt>.
   */
  virtual void get_sub_vector( const Range1D& rng, RTOpPack::MutableSubVector* sub_vec );

  /** \brief Free a mutable explicit view of a sub-vector.
   *
   * The sub-vector view must have been allocated by \c this->get_sub_vector() first.
   *
   * This method has a default implementation which is a companion to the default implementation
   * for <tt>get_sub_vector(...)</tt>.  If <tt>get_sub_vector(...)</tt> is overridden by a subclass then
   * this method must be overridden also!
   *
   *	@param	sub_vec
   *				[in/out] The memory refered to by <tt>sub_vec->values</tt>
   *				and <tt>sub_vec->indices</tt> will be released if it was allocated
   *				and <tt>*sub_vec</tt> will be zeroed out using
   *				<tt>RTOp_mutable_sub_vector_null(sub_vec)</tt>.
   */
  virtual void commit_sub_vector( RTOpPack::MutableSubVector* sub_vec );

  /** \brief Set a specific sub-vector.
   *
   * After this function returns, the corresponding elements in <tt>this</tt> vector object will be
   * set equal to those in the input vector (the post conditions are obvious).
   *
   * Preconditions:<ul>
   * <li> <tt>sub_vec.global_offset + sub_dim <= this->dim()</tt> (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * The default implementation of this operation uses a transformation operator class
   * (see RTOp_TOp_set_sub_vector.h) and calls <tt>apply_op()</tt>.  Be forewarned
   * however, that the operator objects state data (both internal and external) will be 
   * O(<tt>sub_vec.sub_nz</tt>).  For serial applications, this is entirely adequate.  For parallel
   * applications this will be very bad!
   *
   * @param  sub_vec  [in] Represents the elements in the subvector to be set.
   */
  virtual void set_sub_vector( const RTOpPack::SparseSubVector& sub_vec );

  /** \brief Perform a gather or scatter operation with a vector.
   *
     \verbatim

     this = alpha * op(P) * x + beta * this
   \endverbatim
   *
   * The default implementation is based on a transformation or reduction operator
   * (depending if a gather or scatter is being performed).
   */
  virtual void Vp_StMtV(
    value_type                       alpha
    ,const GenPermMatrixSlice        &P
    ,BLAS_Cpp::Transp                P_trans
    ,const Vector                    &x
    ,value_type                      beta
    );

  //@}

  /** @name Overridden from Vector */
  //@{

  /** \brief Default implementation calls <tt>this->sub_view()</tt> (non-<tt>const</tt>) and then
   * performs an cast to <tt>vec_ptr_t</tt>.
   *
   * This function override is actually needed here for another reason.  Without, the
   * override, the non-const version defined in this interface hides the const version
   * defined in Vector.
   */
  vec_ptr_t sub_view( const Range1D& rng ) const;

  //@}

}; // end class VectorMutable

inline
/// <tt>y = alpha * op(P) * x + beta *y</tt>
void Vp_StMtV(
  VectorMutable                    *y	
  ,value_type                      alpha
  ,const GenPermMatrixSlice        &P
  ,BLAS_Cpp::Transp                P_trans
  ,const Vector                    &x
  ,value_type                      beta = 1.0
  )
{
  y->Vp_StMtV(alpha,P,P_trans,x,beta);
}

// ////////////////////////////////////////////////
// Inline members

inline
VectorMutable::vec_mut_ptr_t
VectorMutable::sub_view( const index_type& l, const index_type& u )
{
  return this->sub_view(Range1D(l,u));
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_MUTABLE_HPP
