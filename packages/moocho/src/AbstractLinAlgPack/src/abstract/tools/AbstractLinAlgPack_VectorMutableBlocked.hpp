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

#ifndef VECTOR_WITH_OP_MUTABLE_BLOCK_STD_H
#define VECTOR_WITH_OP_MUTABLE_BLOCK_STD_H

#include <vector>

#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpaceBlocked.hpp"

namespace AbstractLinAlgPack {

/** \brief Concrete subclass for a blocked vector.
 *
 * This is a near optimal implementation for a block vector in many
 * circumstances.
 *
 * A blocked vector is simply the concatonation of two or more vectors to form a larger
 * vector.  Suppose that we are given a list of pointers to vector objects <tt>v[0],v[1],,,v[p]</tt>.
 * Then <tt>this</tt> vector object would be represented as:
 \verbatim

         [ *v[0] ]
         [ *v[1] ]
 this =  [  .    ]
         [ *v[p] ]
 \endverbatim
 * This vector subclass trys to implement all of the vector methods as efficiently as
 * possble.  For example, if a client requests a sub-view that corresponds to a whole constituent
 * sub-vector, then it will return the unadorned vector such as <tt>this->sub_view(1,v[0]->dim()).get() == v[0]</tt>.
 * However, if some sub-view is requested that overlaps two or more constitient vectors, then there is
 * no choice but to return a <tt>%VectorMutableBlocked</tt> object from \c sub_view().
 *
 * There are some situations where this vector subclass will be inefficient.  For example, suppose that
 * the constituent vectors <tt>v[0],v[1],,,v[p]</tt> are parallel vectors with elements owned by unique processes
 * such that an reduction/transformation operator could be applied in parallel to all the constituent vectors in the
 * same time (wall clock) that it takes to apply an operator to one constituent vector.  In this case, this vector
 * subclass could not take advantage of this parallelism and therefore, a more specialized "parallel" block vector
 * class should be created (using threads for instance).  Alternatively, the implementation of this vector
 * subclass could be modified to take advantage of threads if they are available.
 */
class VectorMutableBlocked : virtual public VectorMutable
{
public:

  /** \brief . */
  typedef Teuchos::RCP<const VectorSpaceBlocked>  vec_space_comp_ptr_t;

  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorMutableBlocked(
    VectorMutable::vec_mut_ptr_t*  vecs
    ,const vec_space_comp_ptr_t&         vec_space
    );

  /** \brief Initialize given a list of constituent vectors.
   *
   * Preconditions:<ul>
   * <li> <tt>vec_space.get() != NULL</tt> (throw <tt>std::logic_error</tt>).
   * <li> <tt>vecs[k]->space().is_compatible(*vec_spaces->vector_spaces()[k]) == true</tt>,
   *      for <tt>k = 0...vec_space->num_vector_spaces()</tt>.
   * </ul>
   *
   * Postconditions:<ul>
   * <li> ToDo: Spell out for carefully!
   * </ul>
   *
   * @param  vecs       [in] Array (size <tt>vec_space->num_vector_spaces()</tt>) of smart reference
   *                    counted pointers to the constituent vectors.
   * @param  vec_space  [in] Block vector space object containing constituent vector spaces.
   */
  void initialize(
    VectorMutable::vec_mut_ptr_t*  vecs
    ,const vec_space_comp_ptr_t&         vec_space
    );

  /** \brief Return the underlying block vector space.
   */
  const VectorSpaceBlocked& block_space() const;

  /** \brief Get the kth (zero based) constituent vector.
   */
  const Vector& get_vector(int k) const;

  /** \brief Get the kth (zero based) constituent vector.
   */
  VectorMutable& get_vector(int k);

  //@}

  /** @name Overridden form Vector */
  //@{

  /** \brief . */
  index_type dim() const;
  /** \brief . */
  const VectorSpace& space() const;
  /** \brief . */
  void apply_op(
    const RTOpPack::RTOp& op
    ,const size_t num_vecs, const Vector* vecs[]
    ,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
    ,RTOpPack::ReductTarget *reduct_obj
    ,const index_type first_ele, const index_type sub_dim, const index_type global_offset
    ) const;
  /** \brief . */
  index_type nz() const;
  /** \brief . */
  std::ostream& output(
    std::ostream& out, bool print_dim, bool newline
    ,index_type global_offset 
    ) const;
  /** \brief . */
  value_type get_ele(index_type i) const;
  /** \brief . */
  value_type norm_1() const;
  /** \brief . */
  value_type norm_inf() const;
  /** \brief . */
  value_type inner_product( const Vector& ) const;
  /** \brief . */
  void get_sub_vector( const Range1D& rng, RTOpPack::SubVector* sub_vec ) const;
  /** \brief . */
  void free_sub_vector( RTOpPack::SubVector* sub_vec ) const;

  //@}

  /** @name Overridden from VectorMutable */
  //@{

  /** \brief . */
  vec_mut_ptr_t sub_view( const Range1D& rng );
  /** \brief . */
  void axpy( value_type alpha, const Vector& x );
  /** \brief . */
  VectorMutable& operator=(value_type);
  /** \brief . */
  VectorMutable& operator=(const Vector&);
  /** \brief . */
  void set_ele( index_type i, value_type val );
  /** \brief . */
  void set_sub_vector( const RTOpPack::SparseSubVector& sub_vec );

  //@}

protected:

  /** @name Overridden form Vector */
  //@{
  /** \brief . */
  void has_changed() const;
  //@}

private:

  // ////////////////////////////////////
  // Private types
  
  typedef std::vector<VectorMutable::vec_mut_ptr_t>   vecs_t; ///< Type for list of constituent vectors

  // ////////////////////////////////////
  // Private data members

#ifdef DOXYGEN_COMPILE
  VectorMutable        *vectors;
  VectorSpaceBlocked    *block_vector_space;
#else
  vecs_t                 vecs_;        ///< List of all the vector members.
  vec_space_comp_ptr_t   vec_space_;   ///< overall block vector space.
#endif
  mutable index_type     nz_; ///< > dim() not initalized, < dim() already initialized! 
  mutable value_type     norm_1_, norm_inf_;   ///< < 0, not initialized, > 0 already calculated

  // ////////////////////////////////////
  // Private member functions

  /** \brief . */
  void assert_in_range(int k) const;

  /** \brief . */
  void assert_initialized() const;

  // not defined and not to be called!
  VectorMutableBlocked();
  VectorMutableBlocked(const VectorMutableBlocked&);
  VectorMutableBlocked& operator=(const VectorMutableBlocked&);

}; // end class VectorMutableBlocked

// ////////////////////////////////////
// Inline members

inline
const VectorSpaceBlocked&
VectorMutableBlocked::block_space() const
{
  assert_initialized();
  return *vec_space_;
}

inline
const Vector&
VectorMutableBlocked::get_vector(int k) const
{
  assert_in_range(k);
  return *vecs_[k];
}

inline
VectorMutable&
VectorMutableBlocked::get_vector(int k)
{
  assert_in_range(k);
  return *vecs_[k];
}

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_WITH_OP_MUTABLE_BLOCK_STD_H
