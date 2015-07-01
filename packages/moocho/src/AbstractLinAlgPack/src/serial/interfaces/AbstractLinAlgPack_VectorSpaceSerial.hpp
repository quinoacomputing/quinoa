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

#ifndef VECTOR_SPACE_SERIAL_H
#define VECTOR_SPACE_SERIAL_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Subclass for serial vector space objects that create <tt>VectorMutableDense</tt>
 * vector and <tt>MultiVectorMutableDense</tt> multi-vector objects.
 *
 * The default constructor, copy constructor and assignment operators
 * are allowed since they have the correct semantics.
 */
class VectorSpaceSerial
  : public AbstractLinAlgPack::VectorSpace
{
public:

  /** @name Constructors / initializers */
  //@{

  /** \brief Calls <tt>this->initialize()</tt>.
   */
  VectorSpaceSerial( size_type dim = 0 );

  /** \brief Initialize given the dimension of the vector space.
   *
   * @param  dim   [in] The dimension of the vector space.
   */
  void initialize( size_type dim );

  //@}

  /** @name Overridden from VectorSpece */
  //@{

  /** \brief Returns true if <tt>vec_space.dim() == this->dim()</tt>.
   *
   * The assumption here is that <tt>Vector::get_sub_vector()</tt>
   * and <tt>VectorMutable::get_sub_vector()</tt> can be used to implement
   * all of the methods on an SMP machine in an efficient manner.
   */
   bool is_compatible(const VectorSpace& vec_space) const;
  /// Returns true
  bool is_in_core() const;
  /// Returns 0 if uninitialized
  index_type dim() const;
  /// Returns a <tt>VectorSpaceFactorySerial</tt> object
  space_fcty_ptr_t small_vec_spc_fcty() const;
  /** \brief . */
  space_ptr_t clone() const;
  /// Returns a <tt>VectorMutableDense</tt> object.
  vec_mut_ptr_t create_member() const;
  /// Returns a <tt>MultiVectorMutableDense</tt> object.
  multi_vec_mut_ptr_t create_members(size_type num_vecs) const;
  /** \brief . */
  space_ptr_t sub_space(const Range1D& rng) const;
  /** \brief . */
  space_ptr_t space(
    const GenPermMatrixSlice  &P
    ,BLAS_Cpp::Transp         P_trans
    ) const;
  //@}

private:

  size_type     dim_;

}; // end class VectorSpaceSerial

} // end namespace AbstractLinAlgPack

#endif // VECTOR_SPACE_SERIAL_H
