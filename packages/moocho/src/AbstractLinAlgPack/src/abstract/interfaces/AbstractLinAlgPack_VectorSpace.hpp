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

#ifndef ALAP_VECTOR_SPACE_HPP
#define ALAP_VECTOR_SPACE_HPP

#include "AbstractLinAlgPack_InnerProduct.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface for objects that represent a space for mutable coordinate vectors.
 *
 * This interface acts primarily as an "Abstract Factory" interface for creating \c VectorMutable
 * objects using the \c create_member() method.  A <tt>%VectorSpace</tt> object may also be able
 * to create \c MultiVectorMutable objects which represent a compact collection of vectors.
 * Every application area should be able to define a <tt>%MultiVectorMutable</tt> subclass if
 * it can define a <tt>%VectorMutable</tt> subclass.
 * A secondary role for <tt>%VectorSpace</tt> objects is to test for compatibility of vector spaces
 * (and the vectors and matrix using those spaces) objects using the \c is_compatible() method.
 *
 * Given a <tt>%VectorSpace</tt> object it may also be possible to create sub-spaces using the
 * \c sub_space() method.  This subspace is not a sub-space in the mathematical sense but instead
 * referese to a sub-range of vector elements.
 *
 * Any <tt>%VectorSpace</tt> object can be copied using the \c clone() method.  Therefore,
 * clients have complete control over the lifetime of <tt>%VectorSpace</tt> objects.
 *
 * A <tt>%VectorSpace</tt> object can exist independent from any individual <tt>VectorMutable</tt>
 * (or \c MutiVectorMutable) object; Or, a <tt>%VectorSpace</tt> object can have a lifetime that is
 * dependent on a single <tt>Vector</tt> ( or \c MultiVector) object.  The same interface can
 * serve both roles.
 *
 * Note that a <tt>%VectorSpace</tt> object can create <tt>MultiVectorMutable</tt> objects
 * with any number of column vectors, and <tt>MultiVector::space_rows()</tt> gives a vector space
 * of the dimension.  Therefore, the creation of a multi-vector provides a way for clients to create
 * vector spaces of any arbitrary (although small usually) dimension.  In order to give the client
 * the same ability without having to create a multi-vector, the method <tt>small_vec_spc_fcty()</tt>
 * returns a <tt>VectorSpaceFactory</tt> object that can create vector spaces of small dimension.
 * These created vector spaces can be used by the created <tt>MultiVectorMutable</tt> objects
 * (see <tt>create_members()</tt>).  These vector spaces are also used for the default implementation
 * of <tt>space(P,P_trans)</tt> (see below).  Since a <tt>%VectorSpace</tt> object is not required
 * to create <tt>MultiVectorMutable</tt> objects using <tt>create_members()</tt> a <tt>%VectorSpace</tt>
 * object is also not required to return a <tt>VectorSpaceFactory</tt> object from <tt>small_vec_spc_fcty()</tt>.
 * to be consistent.
 *
 * A vector space is also where the inner product for the space is
 * defined.  It is anticipated that the same implementation of vectors
 * and vector spaces will be reused over and over again with different
 * definitions of the inner product.  Therefore, the inner product is
 * represented as a seperate strategy object.  For example, the same
 * parallel vector implementation can be used with several different
 * inner product definitions.  In some cases, the same inner product
 * stategy object may be able to be used with diffeent vector implemenations
 * (such as the dot product for example).
 *
 * Note that the default copy constructor will transfer the inner product object correctly
 * but the subclass must make sure to copy the inner product object in \c clone() operation.
 * This is little price to pay considering what this design takes care of already for
 * <tt>%VectorSpace</tt> subclasses.
 * However, the default copy constructor will only make a shallow copy of the inner product
 * object but since this object may never be changed, this is perhaps the correct behavior.
 * For any other behavior and the subbclass will have to take care of it.
 *
 * A vector space may also have another advanced feature; it may be able to define other
 * vector spaces based on a gather operation of selected vector elements given a
 * <tt>GenPermMatrixSlice</tt> (<tt>GPMS</tt>) object.  For example, suppose that a
 * <tt>GPMS</tt> object \c P is defined which extracts a (relatively) small number of elements
 * of a vector \c v (from the vector space \c *this) and gathers them into a small vector \c q.
 \verbatim

 q = op(P)*v
 \endverbatim
 * The question is, to what vector space \c Q does the vector \c q belong?
 * The answer is returned by the method <tt>Q = this->space(P,P_trans)</tt>.
 * Theoretically, \c op(P) could be any <tt>GPMS</tt> object and
 * <tt>this->space(P,P_trans)</tt> should always be able
 * to return a non-NULL vector space object.  In reality, the client should only
 * expect <tt>this->space(P,P_trans).get() != NULL</tt> if
 * <tt>q_dim = BLAS_Cpp::rows( P.rows(), P.cols(), P_trans )</tt> is a relatively
 * small number (i.e. less than a few thousand).  If \c q_dim is small, then the vector
 * \c q can always be represented as a local serial vector.  This is not a terribly
 * difficult requirement and any <tt>%VectorSpace</tt> subclass should be able to
 * comply.
 *
 * The returned vector space \c Q must have the property that the scatter operation is
 * also defined.  In other words, it must be that:
 \verbatim

  this->space(P,P_trans)->space(P,trans_not(P_trans))->is_compatible(*this) == true
 \endverbatim
 * This property insures that the gather and scatter operations can be implemented
 * properly.
 */
class VectorSpace
  : public Teuchos::AbstractFactory<VectorMutable>
{
public:

  /// Thrown if vector spaces are incompatible
  class IncompatibleVectorSpaces : public std::logic_error
  {public: IncompatibleVectorSpaces(const std::string& what_arg) : std::logic_error(what_arg) {}};
  /** \brief . */
  typedef Teuchos::RCP<const InnerProduct>          inner_prod_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const VectorSpace>           space_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const VectorSpaceFactory>    space_fcty_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<VectorMutable>               vec_mut_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<MultiVectorMutable>          multi_vec_mut_ptr_t;


  /** @name Constructors / initializers */
  //@{

  /// Calls \c inner_prod()
  VectorSpace( const inner_prod_ptr_t& inner_prod = Teuchos::null );

  /** \brief Initialize with an inner product object.
   *
   * @param  inner_prod  [in] Smart pointer to inner product strategy object.
   *                     If <tt>inner_prod.get()==NULL</tt> then an
   *                     \c InnerProductDot object will be used instead.
   *
   * Postconditions:<ul>
   * <li> [<tt>inner_prod.get() != NULL</tt>] <tt>this->inner_prod().get() == inner_prod.get()</tt>
   * <li> [<tt>inner_prod.get() == NULL</tt>] <tt>dynamic_cast<InnerProductDot*>(this->inner_prod().get()) != NULL</tt>
   * </ul>
   */
  virtual void inner_prod( const inner_prod_ptr_t& inner_prod );

  /** \brief Return the smart pointer to the inner product strategy object.
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * </ul>
   */
  virtual const inner_prod_ptr_t inner_prod() const;

  //@}

  /** @name Pure virtual functions that must be overridden */
  //@{

  /** \brief Create a clone of \c this vector space object.
   *
   * The returned vector space object is expected to be independent from \c this
   * and have a lifetime that extends beyond \c this.  This makes a vector space
   * class a little hander to implement by makes for much better flexibility
   * for the client.  A complete implementation of <tt>%VectorSpace</tt> is not
   * allowed to return \c NULL from this method.
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * </ul>
   */
  virtual space_ptr_t clone() const = 0;

  /** \brief Compare the compatibility of two vector spaces.
   *
   * If this function returns true, then vectors created from
   * either of the vector spaces will be compatible and can
   * be combined in vector operations.
   *
   * Invariants:<ul>
   * <li> [<tt>this->is_compatible(vec_spc) == true</tt>] <tt>vec_spc.is_compatible(*this) == true</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>this->dim() != vec_spc.dim()</tt>] <tt>return == false</tt>
   * </ul>
   */
  virtual bool is_compatible(const VectorSpace& vec_spc ) const = 0;

  /** \brief Return the dimmension of the vector space.
   */
  virtual index_type dim() const = 0;

  /** \brief Create a vector member from the vector space.
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>return->dim() == this->dim()</tt>
   * <li> <tt>return->space().is_compatible(*this) == true</tt>
   * </ul>
   *
   * @return  Returns a smart reference counted pointer to a dynamically
   * allocated vector object.  After construction the values returnd by 
   * <tt>return->get_ele(i)</tt> are unspecified (uninitialized).  This allows for
   * faster execution times.  Note that <tt>&return->space()</tt> does not have to
   * be equal to <tt>this</tt>.
   */
  virtual vec_mut_ptr_t create_member() const = 0;

  //@}

  /** @name Virtual functions with default implementations */
  //@{

  /** \brief Returns true if the vectors are in core.
   *
   * If this function returns true then it means that the vector
   * access functions <tt>Vector::get_sub_vector()</tt> and
   * <tt>VectorMutable::get_sub_vector()</tt> can be safely called and
   * can be expected to be fairly efficient.  If this function does
   * return true then the functions <tt>Vector::get_sub_vector()</tt>,
   * <tt>Vector::free_sub_vector()</tt>,
   * <tt>VectorMutable::get_sub_vector()</tt> and
   * <tt>VectorMutable::commit_sub_vector()</tt> had better be
   * overridden and had better not call
   * <tt>Vector::apply_op(...)</tt>.
   *
   * The default implementation returns <tt>false</tt>
   */
  virtual bool is_in_core() const;

  /** \brief Return a <tt>VectorSpaceFactory</tt> object for the creation of
   * vector spaces with a small dimension.
   *
   * ToDo: Finish documentation!
   *
   * The default implementation returns <tt>return.get() == NULL</tt>.
   *
   */
  virtual space_fcty_ptr_t small_vec_spc_fcty() const;

  /** \brief Create a vector member from the vector space initialized to a scalar.
   *
   * @param  alpha   [in] Scalar that all elements of allocated vector are
   *                 initialized to.
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>return->dim() == this->dim()</tt>
   * <li> <tt>return->space().is_compatible(*this) == true</tt>
   * </ul>
   *
   * @return  Returns a smart reference counted pointer to a dynamically
   * allocated vector object.  After construction the values returnd by 
   * <tt>return->get_ele(i)</tt> are equal to <tt>alpha</tt>.
   * Note that <tt>&return->space()</tt> does not have to
   * be equal to <tt>this</tt>.
   *
   * The default implementation just calls \c create_member() and then
   * assigns <tt>alpha</tt> before returning the smart pointer object.
   */
  virtual vec_mut_ptr_t create_member(const value_type& alpha) const;

  /** \brief Create a set of vector members (a \c MultiVectorMutable) from the vector space.
   *
   * Preconditions:<ul>
   * <li> <tt>num_vecs >= 1</tt> (throw <tt>???</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>return.get() != NULL</tt>] <tt>return->space_cols().is_compatible(*this) == true</tt>
   * <li> [<tt>return.get() != NULL</tt>] <tt>return->space_rows().dim() == num_vecs</tt>
   * <li> [<tt>return.get() != NULL</tt>] <tt>(return->access_by() & MultiVector::COL_ACCESS) == true</tt>
   * </ul>
   *
   * @return  Returns a smart reference counted pointer to a dynamically
   * allocated multi-vector object.  After construction the values returnd by 
   * <tt>return->col(j)->get_ele(i)</tt> are unspecified (uninitialized).  This
   * allows for faster execution times.  Note that <tt>&return->space_cols()</tt>
   * does not have to be equal to <tt>this</tt>.  It is allowed for a vector
   * space implementation to refuse to create multi-vectors and can return
   * \c NULL from this method.
   *
   * The default implementation just returns \c NULL.
   */
  virtual multi_vec_mut_ptr_t create_members(size_type num_vecs) const;

  //@}

  /** \brief Create a transient sub-space of the current vector space.
   *
   * @param  rng  [in] The range of the elements to extract a vector sub-space.
   *
   * Preconditions:<ul>
   * <li> <tt>rng.ubound() <= this->dim()</tt> (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>return.get() != NULL</tt>] <tt>return->dim() == rng->size()</tt>
   * </ul>
   *
   * @return  Returns a smart reference counted pointer to a dynamically
   * allocated vector space object.  Note that the vector object returned
   * by <tt>this->sub_space(rng).create_member()</tt> should be exactly equivalent
   * to the vector returned by
   * <tt>this->create_member()->sub_view(rng)->space()->create_member()</tt>.
   * It is allowed for the implementation to return <tt>return->get() == NULL</tt>
   * for arbitrary values of <tt>rng</tt>.  Only some <tt>rng</tt> ranges may be allowed
   * but they will be appropriate for the application at hand.  However, a
   * very good implementation should be able to accomidate any valid <tt>rng</tt>
   * that meets the basic preconditions.
   *
   * Note that if two vector space objects <tt>X</tt> and <tt>Y</tt> are compatible
   * (i.e. <tt>X.is_compatible(Y) == true</tt>, then it is also expected that
   * <tt>X.sub_space(rng)->is_compatible(*Y.sub_space(rng))</tt> will also be \c true.
   * However, in general it can not be expected that
   * <tt>X.sub_space(rng1)->is_compatible(*X.sub_space(rng2)) == true</tt>, even if
   * <tt>rng1.size() == rng2.size()</tt>.  For serial vectors, it may
   * be but for parallel vectors it will most certainly not be.  Therefore, in
   * general, don't assume that arbitrary subsets of the vector spaces will be
   * compatible, even if the sizes of these subspaces are the same.
   *
   * The default implementation uses the subclass \c VectorSpaceSubSpace to
   * represent any arbitrary sub-space but this can be very inefficient if the
   * sub-space is very small compared this this full vector space.
   */
  virtual space_ptr_t sub_space(const Range1D& rng) const;

  /// Inlined to call <tt>this->sub_space(Range1D(il,iu))</tt>.
  space_ptr_t sub_space( const index_type il, const index_type iu ) const;

  /** \brief Create a vector space for vector to gather the elements into.
   *
   * @param  P        [in] A <tt>GenPermMatrixSlice</tt> object specifying the map.
   * @param  P_trans  [in] Determines if <tt>P</tt> is transposed or not.
   *
   * Preconditions:<ul>
   * <li> ???
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>return.get()!=NULL</tt>]
   *   return->space(P,trans_not(P_trans))->is_compatible(*this) == true
   * </ul>
   *
   * @return  Returns a smart reference counted pointer to a dynamically
   * allocated vector space object.  
   *
   * ToDo: Finish documentation!
   *
   * The default implementation returns <tt>return.get() == NULL</tt>.
   */
  virtual space_ptr_t space(
    const GenPermMatrixSlice  &P
    ,BLAS_Cpp::Transp         P_trans
    ) const;

  //@}
  
  /** @name Overridden from AbstractFactory */
  //@{

  /// Just calls <tt>this->create_member()</tt> by default!
  obj_ptr_t create() const;

  //@}

private:
#ifdef DOXYGEN_COMPILE
  const InnerProduct   *inner_prod;
#else
  inner_prod_ptr_t     inner_prod_;
#endif

}; // end class VectorSpace

// ////////////////////////////////////////////////
// Inline members

inline
VectorSpace::space_ptr_t
VectorSpace::sub_space( const index_type il, const index_type iu ) const
{
  return this->sub_space(Range1D(il,iu));
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_SPACE_HPP
