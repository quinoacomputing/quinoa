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

#ifndef ALAP_VECTOR_HPP
#define ALAP_VECTOR_HPP

#include <iosfwd>

#include "AbstractLinAlgPack_Types.hpp"
#include "RTOpPack_RTOpT.hpp"

namespace AbstractLinAlgPack {

/** \brief Apply a reduction/transformation,operation over a set of vectors:
 * <tt>op(op(v[0]...v[nv-1],z[0]...z[nz-1]),(*reduct_obj)) -> z[0]...z[nz-1],(*reduct_obj)</tt>.
 *
 * The logical vector passed to the
 * <tt>op\ref RTOpPack::RTOp::apply_op ".apply_op(...)"</tt>
 * method is: \verbatim

 v(k + global_offset) = this->get_ele(first_ele + k - 1)
 , for k = 1 ... sub_dim
 \endverbatim
 *
 * where <tt>v</tt> represents any one of the input or input/output
 * vectors.  The situation where <tt>first_ele == 1</tt> and
 * <tt>global_offset > 1</tt> corresponds to the case where the
 * vectors represent consitituent vectors in a larger aggregrate
 * vector.  The situation where <tt>first_ele > 1</tt> and
 * <tt>global_offset == 0</tt> is for when a sub-view of the vectors
 * are being treated as full vectors.  Other combinations of these
 * arguments are also possible.
 *
 * Preconditions:<ul>
 * <li> [<tt>num_vecs > 0</tt>] <tt>vecs[k]->space().isCompatible(this->space()) == true</tt>
 *          , for <tt>k = 0...num_vecs-1</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
 * <li> [<tt>num_targ_vecs > 0</tt>] <tt>targ_vecs[k]->space().isCompatible(this->space()) == true</tt>
 *          , for <tt>k = 0...num_targ_vecs-1</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
 * <li> <tt>1 <= first_ele <= this->dim()</tt> (throw <tt>std::out_of_range</tt>)
 * <li> <tt>global_offset >= 0</tt> (throw <tt>std::invalid_argument</tt>)
 * <li> <tt>sub_dim - (first_ele - 1) <= this->dim()</tt> (throw <tt>std::length_error</tt>).
 * </ul>
 *
 * @param  op	[in] Reduction/transformation operator to apply over each sub-vector
 *				and assemble the intermediate targets into <tt>reduct_obj</tt> (if
 *              <tt>reduct_obj != NULL</tt>).
 * @param  num_vecs
 *				[in] Number of nonmutable vectors in <tt>vecs[]</tt>.
 *              If <tt>vecs==NULL</tt> then this argument is ignored but should be set to zero.
 * @param  vecs
 *				[in] Array (length <tt>num_vecs</tt>) of a set of pointers to
 *				nonmutable vectors to include in the operation.
 *				The order of these vectors is significant to <tt>op</tt>.
 * @param  num_targ_vecs
 *				[in] Number of mutable vectors in <tt>targ_vecs[]</tt>.
 *              If <tt>targ_vecs==NULL</tt>	then this argument is ignored but should be set to zero.
 * @param  targ_vecs
 *				[in] Array (length <tt>num_targ_vecs</tt>) of a set of pointers to
 *				mutable vectors to include in the operation.
 *				The order of these vectors is significant to <tt>op</tt>.
 *				If <tt>targ_vecs==NULL</tt> then <tt>op</tt> is called with no mutable vectors.
 * @param  reduct_obj
 *				[in/out] Target object of the reduction operation.
 *				This object must have been created by the <tt>op.reduct_obj_create_raw(&reduct_obj)</tt>
 *              function first.  The reduction operation will be added to <tt>(*reduct_obj)</tt> if
 *              <tt>(*reduct_obj)</tt> has already been through a reduction.  By allowing the info in
 *              <tt>(*reduct_obj)</tt> to be added to the reduction over all of these vectors, the reduction
 *              operation can be accumulated over a set of abstract vectors	which can be useful for implementing
 *              composite vectors for instance.  If <tt>op.get_reduct_type_num_entries(...)</tt> returns
 *              <tt>num_values == 0</tt>, <tt>num_indexes == 0</tt> and <tt>num_chars == 0</tt> then
 *              <tt>reduct_obj</tt> should be set to <tt>NULL</tt> and no reduction will be performed.
 * @param  first_ele
 *				[in] (default = 1) The index of the first element in <tt>this</tt> to be included.
 * @param  sub_dim
 *              [in] (default = 0) The number of elements in these vectors to include in the reduction/transformation
 *              operation.  The value of <tt>sub_dim == 0</tt> means to include all available elements.
 * @param  global_offset
 *				[in] (default = 0) The offset applied to the included vector elements.
 *
 */
void apply_op(
  const RTOpPack::RTOp       &op
  ,const size_t              num_vecs
  ,const Vector*             vecs[]
  ,const size_t              num_targ_vecs
  ,VectorMutable*            targ_vecs[]
  ,RTOpPack::ReductTarget    *reduct_obj
  ,const index_type          first_ele      = 1
  ,const index_type          sub_dim        = 0
  ,const index_type          global_offset  = 0
  );

/** \brief Abstract interface for immutable, finite dimensional, coordinate vectors {abstract}.
 *
 * This interface contains a mimimal set of operations.  The main feature
 * of this interface is the operation <tt>apply_op()</tt>.
 * Almost every standard (i.e. BLAS) and non-standard element-wise operation that
 * can be performed on a set of coordinate vectors without changing (mutating)
 * the vectors can be performed through reduction operators.  More standard
 * vector operations could be included in this interface and allow
 * for specialized implementations but in general, assuming the
 * sub-vectors are large enough, such implementations
 * would not be significantly faster than those implemented through
 * reduction/transformation operators.  There are some operations however
 * that can not always be efficiently with reduction/transforamtion operators
 * and a few of these important methods are included in this interface.  The
 * <tt>apply_op()</tt> method allows to client to specify a sub-set
 * of the vector elements to include in reduction/transformation operation.
 * This greatly increases the generality of this vector interface as vector
 * objects can be used as sub objects in larger composite vectors and sub-views
 * of a vector can be created.
 *
 * This interface allows clients to create sub-views of a vector using <tt>sub_view()</tt>
 * that in turn are fully functional <tt>%Vector</tt> objects.  This functionality
 * is supported by default by using a default vector subclass <tt>VectorSubView</tt> which
 * in turn calls <tt>apply_op()</tt> but the client need not ever worry about
 * how this is done.
 *
 * This interface also allows a client to extract a sub-set of elements in an
 * explicit form as an <tt>RTOpPack::SubVector object using the method \c get_sub_vector().
 * In general, this is very bad thing to do and should be avoided at all costs.
 * However, there are some applications where this is needed and therefore it is
 * supported.  The default implementation of this method uses a reduction/transformation
 * operator with <tt>apply_op()</tt> in order to extract the needed elements.
 *
 * In order to create a concreate subclass of this interface, only two
 * methods must be overridden.  The <tt>space()</tt> method must be overridden which in turn
 * requires defining a concreate <tt>VectorSpace</tt> class (which has only two pure virtual
 * methods).  And, as mentioned above, the <tt>apply_op()</tt> method must be overridden
 * as well.
 *
 * The fact that this interface defines <tt>space()</tt> which returns a <tt>VectorSpace</tt> object
 * (which in turn can create mutable vectors) implies that for every possible vector object,
 * it is possible to associate with it a mutable vector object that can be the target
 * of transformation operations.  This is not a serious limitation.  For any
 * application area, mutable vectors should be able to defined and should be
 * usable with the non-mutable vectors.
 *
 * This interface includes methods for the common vector norms: \c norm_1(),
 * \c norm_2(), \c norm_inf().  The default implementation of this class uses reduction
 * operator classes (See RTOp_ROp_norms.h) and caches away the values of the
 * norms that are computed since it is common that the norms will be accessed many
 * times before a vector is changed.  The operations in any subclass that modifies
 * the underlying vector must call the method <tt>this-></tt>has_changed() in order
 * to alert this implementation that the norms are no longer valid.
 *
 * ToDo: Add example code!
 */
class Vector {
public:

  /** \brief . */
  typedef Teuchos::RCP<const Vector>   vec_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<VectorMutable>  vec_mut_ptr_t;

  /** \brief . */
  friend
  void apply_op(
    const RTOpPack::RTOp       &op
    ,const size_t              num_vecs
    ,const Vector*             vecs[]
    ,const size_t              num_targ_vecs
    ,VectorMutable*            targ_vecs[]
    ,RTOpPack::ReductTarget    *reduct_obj
    ,const index_type          first_ele
    ,const index_type          sub_dim
    ,const index_type          global_offset
    );

  /** \brief . */
  Vector();
  /** \brief . */
  virtual ~Vector() {}

  /** @name Pure virtual methods (must be overridden by subclass) */
  //@{

  /** \brief Return the vector space that this vector belongs to.
   *
   * Note that the vectors space object returned is specifically bound to this
   * vector object.  The vector space object returned should only be considered
   * to be transient and may become invalid if <tt>this</tt> is modified in some way
   * (though some other interface).
   */
  virtual const VectorSpace& space() const = 0;

protected:

  /** \brief Apply a reduction/transformation,operation over a set of vectors:
   * <tt>op(op(v[0]...v[nv-1],z[0]...z[nz-1]),(*reduct_obj)) -> z[0]...z[nz-1],(*reduct_obj)</tt>.
   *
   * The vector <tt>this</tt> that this method is called on is
   * assumed to be one of the vectors in
   * <tt>v[0]...v[nv-1],z[0]...z[nz-1]</tt>.  This method is not
   * called directly by the client but is instead
   * <tt>TSFCore::applyOp()</tt>.
   *
   * See the documentation for the method <tt>AbstractLinAlgPack::apply_op()</tt>
   * for a description of what this method does.
   */
  virtual void apply_op(
    const RTOpPack::RTOp       &op
    ,const size_t              num_vecs
    ,const Vector*             vecs[]
    ,const size_t              num_targ_vecs
    ,VectorMutable*            targ_vecs[]
    ,RTOpPack::ReductTarget    *reduct_obj
    ,const index_type          first_ele
    ,const index_type          sub_dim
    ,const index_type          global_offset
    ) const = 0;

  //@}

public:

  /** @name Miscellaneous virtual methods with default implementations */
  //@{

  /** \brief Return the dimension of this vector.
   *
   * It is allowed for a vector to return a dimension of <tt>0</tt> in which case
   * the vector should be considered uninitialized in which the client should
   * not call any of the member functions (except space()).  The default implementation
   * returns <tt>this->space().dim()</tt>.
   */
  virtual index_type dim() const;

  /** \brief Return the number of nonzero elements in the vector.
   *
   * The default implementation just uses a reduction operator
   * with the <tt>apply_op()</tt> method (See
   * RTOp_ROp_num_nonzeros.h).
   */
  virtual index_type nz() const;

  /** \brief Virtual output function.
    *
    * The default implementation just uses get_sub_vector<tt>(...)</tt> to convert to
    * a dense vector and then prints this.
    *
    * ToDo: Finish documentation!
    *
    * @param  out        [in/out] Receives the output.
    * @param  print_dim  [in] (default = true) If true, then the dimension is printed
    *                    first followed by a newline.
    * @param  newline    [in] (default = true) If true, then a newline is printed after
    *                    the last element is printed.
    * @param  global_offset
    *                    [in] (default = 0) The offset added to the vector element
    *                    indexes when they are printed.
    */
  virtual std::ostream& output(
    std::ostream& out, bool print_dim = true, bool newline = true
    ,index_type global_offset = 0
    ) const;

  /** \brief Create a clone of this vector objet.
   *
   * The vector object returned in a smart reference counted pointer to a functional copy of
   * the current vector object.  The vector object <tt>this</tt> and the vector returned by
   * this method can be modified independently.
   *
   * The default implementation of this function calls on <tt>this->space().create_member()</tt> and
   * then copies over the elements from <tt>this</tt> using <tt>operator=()</tt>.
   */
  virtual vec_mut_ptr_t clone() const;

  /** \brief Fetch an element in the vector.
   *
   * Preconditions:<ul>
   * <li> <tt>1 <= i <= this->dim()</tt> (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * The default implementation uses a C reduction operator class
   * (See RTOp_ROp_get_ele.h C).
   *
   * @param  i  [in]  Index of the element value to get.
   */
  virtual value_type get_ele(index_type i) const;

  /** \brief Create an abstract view of a vector object .
   *
   * This is only a transient view of a sub-vector that is to be immediately used
   * and then released by <tt>RCP<></tt>.
   *
   * It is important to understand what the minimum postconditions are for the sub vector objects
   * returned from this method.  If two vector objects <tt>x</tt> and <tt>y</tt> are compatible (possibly of
   * different types) it is assumed that <tt>*x.sub_view(rng)</tt> and <tt>*y.sub_view(rng)</tt>
   * will also be compatible vector objects no mater what range <tt>rng</tt> represents.  However,
   * if <tt>i1 < i2 < i3 < i4</tt> with <tt>i2-i1 == i4-i3</tt>, then in general, one can not expect
   * that the vector objects <tt>*x.sub_view(Range1D(i2,i1))</tt> and
   * <tt>*x.sub_view(Range1D(i4,i5))</tt> will be compatible objects.  For some vector
   * implementaions they may be (i.e. serial vectors) but for others they most
   * certainly will not be (i.e. parallel vectors).  This limitation must be kept in
   * mind by all vector subclass implementors and vector interface clients.
   *
   * Preconditions:<ul>
   * <li> [<tt>!rng.full_range()</tt>] <tt>(rng.ubound() <= this->dim()) == true</tt>
   *      (<tt>throw std::out_of_range</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>[return.get() != NULL] return->get_ele(i) == this->get_ele(i+rng.lbound()-1)</tt>
   *       , for <tt>i = 1...rng.size()</tt>.
   * </ul>
   *
   * @param  rng  [in] The range of the elements to extract the sub-vector view.  It
   *              is allowed for <tt>rng.full_range() == true</tt> in which case it implicitly
   *              treated as <tt>rng = [1,this->dim()]</tt>.
   * 
   * @return  Returns a smart reference counted pointer to a view of the requested
   * vector elements.  It is allowed for subclasses to return  <tt>return->get() == NULL</tt>
   * for some selections
   * of <tt>rng</tt>.  Only some <tt>rng</tt> ranges may be allowed but they will be appropriate for the
   * application at hand.  However, a very good implementation should be able to
   * accommodate any valid <tt>rng</tt> that meets the basic preconditions.  The default
   * implementation uses the subclass \c VectorSubView to represent any arbitrary
   * sub-view but this can be inefficient if the sub-view is very small compared this this
   * full vector space but not necessarily.
   */
  virtual vec_ptr_t sub_view( const Range1D& rng ) const;

  //@}

  /** \brief Inline member function that simply calls <tt>this->sub_view(Range1D(l,u))</tt>.
   */
  vec_ptr_t sub_view( const index_type& l, const index_type& u ) const;

  /** @name Vector norms */
  //@{

  /** \brief One norm. <tt>||v||_1 = sum( |v(i)|, i = 1,,,this->dim() )</tt>
   */
  virtual value_type norm_1() const;
  /** \brief Two norm. <tt>||v||_2 = sqrt( sum( v(i)^2, i = 1,,,this->dim() ) )</tt>
   */
  virtual value_type norm_2() const;
  /** \brief Infinity norm.  <tt>||v||_inf = max( |v(i)|, i = 1,,,this->dim() )</tt>
   */
  virtual value_type norm_inf() const;
  
  //@}

  /** @name Inner product */
  //@{

  /** \brief Return the inner product of <tt>*this</tt> with <tt>v</tt>.
   *
   * @return Returns <tt>this->space().inner_prod()->inner_prod(*this,v)</tt>
   */
  virtual value_type inner_product( const Vector& v ) const;

  //@}

  /** @name Explicit sub-vector access */
  //@{

  /** \brief Get a non-mutable explicit view of a sub-vector.
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
   * If a subclass does override this method, it must also override <tt>release_sub_vector(...)</tt>
   * which has a default implementation which is a companion to this method's default
   * implementation.
   *
   * @param  rng      [in] The range of the elements to extract the sub-vector view.
   * @param  sub_vec  [in/out] View of the sub-vector.  Prior to the
   *                  first call <tt>RTOp_sub_vector_null(sub_vec)</tt> must
   *                  have been called for the correct behavior.  Technically
   *                  <tt>*sub_vec</tt> owns the memory but this memory can be freed
   *                  only by calling <tt>this->free_sub_vector(sub_vec)</tt>.
   */
  virtual void get_sub_vector( const Range1D& rng, RTOpPack::SubVector* sub_vec ) const;

  /** \brief Free an explicit view of a sub-vector.
   *
   * The sub-vector view must have been allocated by this->get_sub_vector() first.
   *
   * This method has a default implementation which is a companion to the default implementation
   * for <tt>get_sub_vector(...)</tt>.  If <tt>get_sub_vector(...)</tt> is overridden by a subclass then
   * this method must be overridden also!
   *
   *	@param	sub_vec
   *				[in/out] The memory refered to by <tt>sub_vec->values</tt>
   *				and <tt>sub_vec->indices</tt> will be released if it was allocated
   *				and <tt>*sub_vec</tt> will be zeroed out using
   *				<tt>RTOp_sub_vector_null(sub_vec)</tt>.
   */
  virtual void free_sub_vector( RTOpPack::SubVector* sub_vec ) const;

  //@}

  /** \brief Must be called by any vector subclass that modifies this vector
   * object!
   *
   * The way to use this method by subclasses is to call it when ever
   * there is a chance that the vector may have changed.  Therefore, this
   * method should be called in every non-const member function in every
   * subclass.  This is a little bit of a pain but overall this is a good
   * design in that it allows for efficient cacheing of information for
   * multiple retreval.  For example, if the subclass <tt>SomeVector</tt> has cashed
   * data and has a method <tt>SomeVector::foo()</tt> may modify the
   * vector then <tt>SomeVector</tt> should override the method <tt>has_changed()</tt> and its
   * implementation should look someting likde like this!
   \verbatim
   void SomeVector::has_changed()
   {
       BaseClass::has_changed(); // Called on most direct subclass that
                               // has overridden this method as well.
      ...  // Reinitialize your own cached data to uninitialized!
   }
   \endverbatim
   */
  virtual void has_changed() const;

protected:

  /** @name Protected helper functions */
  //@{

  /** \brief This method usually needs to be called by subclasses at the
   * end of the <tt>apply_op()</tt> method implementation to
   * insure that <tt>has_changed()</tt> is called on the transformed
   * vector objects.
   */
  virtual void finalize_apply_op(
    const size_t num_targ_vecs, VectorMutable** targ_vecs
    ) const;

  //@}

private:

  mutable index_type  num_nonzeros_;  ///< < 0 == not initialized, > 0 == already calculated
  mutable value_type  norm_1_, norm_2_, norm_inf_;  ///< < 0 == not initialized, > 0 == already calculated

}; // end class Vector

// ////////////////////////////////////////////////
// Inline functions

inline
Vector::vec_ptr_t
Vector::sub_view( const index_type& l, const index_type& u ) const
{
  return this->sub_view(Range1D(l,u));
}

} // end namespace AbstractLinAlgPack

#endif  // ALAP_VECTOR_HPP
