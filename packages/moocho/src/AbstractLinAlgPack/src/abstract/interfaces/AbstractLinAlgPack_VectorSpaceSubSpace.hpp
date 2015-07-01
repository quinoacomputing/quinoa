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

#ifndef VECTOR_SPACE_SUB_SPACE_H
#define VECTOR_SPACE_SUB_SPACE_H

#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Concrete subclass for a default sub-space of a vector.
 *
 * There is not much to this subclass.  It basically implements all
 * of its methods based on the external VectorSpace interface to
 * implement is_compatible() and sub_space() and and relys
 * on a default subclass VectorMutableSubView to implement
 * create_member(). 
 *
 * The default constructor, copy constructor and assignment operator
 * functions are allowed and have the correct behavior.
 */
class VectorSpaceSubSpace : public virtual VectorSpace {
public:

  /** \brief Constructs to uninitialized.
   *
   * Postconditions: see \c set_uninitialized().
   */
  VectorSpaceSubSpace();

  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorSpaceSubSpace( const space_ptr_t& full_space, const Range1D& rng );

  /** \brief Initialize.
   *
   * Constructs a sub-space of the vector space this = space.sub_space(rng).
   *
   * Preconditions:<ul>
   * <li> <tt>full_space.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>).
   * <li> [rng.full_range() == false</tt>] <tt>rng.lbound() <= full_space->dim()</tt> (throw <tt>std::out_of_range</tt>).
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [&& rng.full_range() == true</tt>] <tt>this->dim() == full_space->dim()</tt>
   * <li> [&& rng.full_range() == false</tt>] <tt>this->dim() == rng.size()</tt>
   * </ul>
   *
   * @param  full_space  [in] The original full vector space (must be <tt>full_space.get() != NULL</tt>).
   * @param  rng         [in] The range of element that <tt>this</tt> vector sub-space will represent.
   */
  void initialize( const space_ptr_t& full_space, const Range1D& rng );

  /** \brief Set uninitilized.
   *
   * Postconditions:<ul>
   * <li> <tt>this->dim() == 0</tt>
   * <li> <tt>this->create_member().get() == NULL</tt>
   * </ul>
   */
  void set_uninitialized();

  /** \brief . */
  const space_ptr_t& full_space() const;

  /** \brief . */
  const Range1D& rng() const;

  /// Validate rng
  void validate_range( const Range1D& rng ) const;

  /** @name Overridden from VectorSpace */
  //@{

  /** \brief . */
  bool is_compatible(const VectorSpace& ) const;
  /** \brief . */
  bool is_in_core() const;
  /** \brief . */
  index_type dim() const;
  /** \brief . */
  vec_mut_ptr_t create_member() const;
  /** \brief . */
  space_ptr_t clone() const;
  /** \brief . */
  space_ptr_t sub_space(const Range1D& rng) const;

  //@}

private:

  space_ptr_t     full_space_;   ///< If space_.get() == NULL, then uninitalized (dim == 0)
  Range1D         rng_;          ///< The range of elements from this space to represent!

}; // end class VectorSpaceSubSpace

// //////////////////////////////
// Inline members

inline
VectorSpaceSubSpace::VectorSpaceSubSpace()
  : rng_(Range1D::Invalid)
{}

inline
const VectorSpace::space_ptr_t& VectorSpaceSubSpace::full_space() const
{
  return full_space_;
}

inline
const Range1D& VectorSpaceSubSpace::rng() const
{
  return rng_;
}

#ifndef TEUCHOS_DEBUG
inline
void VectorSpaceSubSpace::validate_range(const Range1D&) const
{}
#endif

} // end namespace AbstractLinAlgPack

#endif // VECTOR_SPACE_SUB_SPACE_H
