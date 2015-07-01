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

#ifndef VECTOR_SPACE_COMPOSITE_STE_H
#define VECTOR_SPACE_COMPOSITE_STE_H

#include <vector>

#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief \c VectorSpace subclass for the composite of one or more <tt>%VectorSpace</tt>
 * objects.
 *
 * This subclass allows <tt>%VectorSpace</tt> objects to be built out of one or more
 * other vector space objects.  This is essentially the product of vector
 * spaces.  The specific type of vector created by a \c VectorSpaceBlocked
 * object is of type \c VectorMutableBlocked but the client need not
 * ever know this or deal with this type directly.
 *
 * For example, suppose you have \c p vector spaces <tt>V[k]</tt> for <tt>k = 0...p-1</tt>
 * and want to form a concatenated vector space \c Z containing all of these vector spaces
 * stacked on top of each other to form:
 \verbatim

     [ V[0]   ]
 Z = [ V[1]   ]
     [ .      ]
     [ V[p-1] ]
 \endverbatim
 * Such a vector space can be constructed as shown if the following function:
 \code
 void f( const VectorSpace::space_ptr_t V[], int p, VectorSpaceBlocked* Z )
 {
     Z->initialize( V, p );
 }
 \endcode
 * Once a <tt>%VectorSpaceBlocked</tt> object is initialized, it can be used just like
 * any other <tt>%VectorSpace</tt> object.  The method \c create_member() will create
 * \c VectorWithOpMutableCompoisteStd objects containing members from the constitient
 * vector spaces.
 *
 * There are several methods that can be used by clients that need to work with the
 * individual constituent vector spaces.  The method \c num_vector_spaces() give the
 * number of constitient vector spaces while \c vector_spaces() returns an array
 * of the constitient vector spaces passed to \c initialize().  Some other useful
 * utility methods are also defined.  The method \c vector_spaces_offsets() returns
 * an array that gives the offset of each constitient vector in the overall composite
 * vector.  For example, the <tt>ith</tt> (zero based) vector space <tt>this->%vector_spaces()[i]</tt>
 * owns the element indexes <tt>this->%vector_spaces_offsets()[i]+1</tt> to
 * <tt>this->%vector_spaces_offsets()[i+1]</tt>.  Determining which constitient vector
 * space owns a element index can be determined by calling \c get_vector_space_position().
 *
 * The default assignment operator is allowed since it has the correct semantics.
 * The default copy constructor is also allowed but only performs a shallow copy
 * of the constituent vector space objects.  If you want to copy the constituent
 * vector space objects also you need to use the \c clone() method.
 * The default constructor is not allowed (declared private) to avoid accidents.
 */
class VectorSpaceBlocked : public VectorSpace {
public:

  /** \brief Construct the vector space object.
   *
   * Calls <tt>this->initialize()</tt>.
   */
  VectorSpaceBlocked(
    const VectorSpace::space_ptr_t          vector_spaces[]
    ,int                                    num_vector_spaces
    ,const VectorSpace::space_fcty_ptr_t   &small_vec_spc_fcty = Teuchos::null
    );
  
  /** \brief Initialize with a set of vector space objects.
   *
   * Postconditions:<ul>
   * <li> <tt>this->num_vector_spaces() == num_vector_spaces</tt>
   * <li> <tt>this->vector_spaces()[i].get() == vector_spaces[i].get(), i = 0...num_vector_spaces-1</tt>
   * <li> <tt>this->dim() == sum{ vector_spaces[i]->dim(), i = 0...num_vector_spaces-1 }</tt>
   * <li> <tt>this->sub_space(rng[i]).get() == vector_spaces[i].get()</tt>,
     *      <tt>this->vector_spaces_offsets()[i] = rng[i].lbound()-1</tt>,
   *      for <tt>i = 0...num_vector_spaces-1</tt>, <tt>rng[0] = Range1D(1,vector_spaces[0]->dim());</tt>
   *      <tt>rng[i] = Range1D(1,vector_spaces[i]->dim()) + rng[i-1].ubound(), i = 0...num_vector_spaces-1</tt>
     * <li> <tt>this->dim() == this->vector_spaces_offsets()[num_vector_spaces]</tt>
   * <li> <tt>this->small_vec_spc_fcty().get() == small_vec_spc_fcty.get()</tt>
   * </ul>
    *
   * @param  vector_spaces
   *               [in] Array (length \c num_vector_spaces) that give the
   *                vector space objects that compose this composite vector space.
   *                The quantity <tt>vector_spaces[i].get()</tt>, for <tt>i = 0...num_vector_spaces-1</tt>
   *                gives the <tt>ith</tt> vector_space object.
   *                <tt>vector_spaces</tt> an be NULL if <tt>num_vector_spaces == 0</tt>.
   * @param  num_vector_spaces
   *               [in] Number of constitient vector spaces.  Can be zero (empty vector space).
   * @param  small_vec_spc_fcty
   *               [in] Defines the vector space factory returned by <tt>this->small_vec_spc_fcty()</tt>
   */
  void initialize(
    const VectorSpace::space_ptr_t         vec_spaces[]
    ,int                                   num_vector_spaces
    ,const VectorSpace::space_fcty_ptr_t   &small_vec_spc_fcty = Teuchos::null
    );

  /// Return the value of <tt>num_vec_spaces</tt> passed to <tt>this->initialize()</tt>.
  int num_vector_spaces() const;

  /** \brief . */
    /* Returns an array (length <tt>this->num_vector_spaces()</tt>) consistent with
     * <tt>vec_spaces</tt> passed into <tt>this->initialize()</tt>.
     */
  const VectorSpace::space_ptr_t*  vector_spaces() const;

  /** \brief Returns array (length <tt>this->num_vector_spaces() + 1</tt>) of offsets for each constitient
     * vector space in the composite vector.
     */
  const index_type* vector_spaces_offsets() const;

  /** \brief Get the position of the vector space object and its offset into the composite
   * vector of the vector space object that owns the <tt>ith</tt> index in the composite
   * vector.
   *
   * Preconditions:<ul>
   * <li> <tt>1 <= i <= this->dim()</tt> (throw <tt>std::out_of_range</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>kth_global_offset + 1 <= i <= kth_global_offset + this->vector_spaces()[kth_vector_space]->dim()</tt>
   * </ul>
   *
   * @param  i    [in] The index of the element to find the vector space object for.
   * @param  kth_vector_space
   *              [out] The index for <tt>this->vector_spaces()[kth_vector_space]</tt> that owns the element <tt>i</tt>.
   * @param  kth_global_offset
   *              [out] The global offset for <tt>this->vector_spaces()[kth_vector_space]</tt> in the composite
   *              vector.
   */
  void get_vector_space_position( index_type i, int* kth_vector_space, index_type* kth_global_offset ) const;

  /** @name Overridden from VectorSpace */
  //@{
  
  /** \brief Returns \c true if same type and has compatible vector spaces.
     *
     * This method will only return true if all of the following are true:
     * <ul>
     * <li> <tt>(vs = dynamic_cast<const VectorSpaceCompoisteStd*>(&vec_space)) != NULL</tt> </li>
     * <li> <tt>this->num_num_vector_spaces() == vs->num_vector_spaces()</tt> </li>
     * <li> <tt>this->vector_spaces()[k]->is_compatible(*vs->vector_spaces()[k]<tt>
     *      , for <tt>k = 0 ... this->num_vector_spaces()-1</tt> </li>
     * </ul>
     * Otherwise, this method will return false.
   */
   bool is_compatible(const VectorSpace& vec_space) const;
  /** \brief . */
  index_type dim() const;
  /** \brief . */
  space_fcty_ptr_t small_vec_spc_fcty() const;
  /** \brief . */
  vec_mut_ptr_t create_member() const;
  /** \brief . */
  multi_vec_mut_ptr_t create_members(size_type num_vecs) const;
  /** \brief . */
  space_ptr_t clone() const;
  /** \brief . */
  space_ptr_t sub_space(const Range1D& rng) const;

  //@}

private:

  // ////////////////////////////////////////////
  // Private types

  typedef std::vector<VectorSpace::space_ptr_t>   vector_spaces_t;
  typedef std::vector<index_type>                 vec_spaces_offsets_t;

  // ////////////////////////////////////////////
  // Private data members

#ifdef DOXYGEN_COMPILE
  VectorSpace           *vector_spaces;
  VectorSpaceFactory    *small_vec_spc_fcty;
#else
  vector_spaces_t       vector_spaces_;
  // Pointer to the list of vector spaces.

  vec_spaces_offsets_t  vec_spaces_offsets_;
  // vec_spaces_offset_[k] gives the offset for vector_spaces_[k] in the 
  // global composite vector.  vec_spaces_offset_[vector_spaces_.size()] gives
  // the total dimension of the global composite vector.

  VectorSpace::space_fcty_ptr_t   small_vec_spc_fcty_;
  // Vector space factory for small_vec_spc_fcty
#endif

  // Not defined and not to be called
  VectorSpaceBlocked();

}; // end class VectorSpaceBlocked

// /////////////////////////////////////////////////////////
// Inline members

inline
int VectorSpaceBlocked::num_vector_spaces() const
{
  return vector_spaces_.size();
}

inline
const VectorSpace::space_ptr_t*
VectorSpaceBlocked::vector_spaces() const
{
  return &vector_spaces_[0];
}

inline
const index_type* VectorSpaceBlocked::vector_spaces_offsets() const
{
  return &vec_spaces_offsets_[0];
}

} // end namespace AbstractLinAlgPack

#endif // VECTOR_SPACE_COMPOSITE_STE_H
