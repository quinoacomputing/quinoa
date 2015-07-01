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

#ifndef ALAP_PERMUTATION_H
#define ALAP_PERMUTATION_H

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract interface to permutation matrices.
 *
 * A \c Permutation object is used to permute the elements within
 * a vector.  It is not a general linear operator since it does not
 * map between vector spaces.  It only permutes elements within the same
 * vector space.
 */
class Permutation {
public:

  /** \brief . */
  virtual ~Permutation() {}

  /** @name Vector space */
  //@{
  
  /** \brief Return a reference to a vector space object that this permutation is compatible with.
   */
  virtual const VectorSpace& space() const = 0;
  
  //@}

  /** @name Information */
  //@{

  /// Return the dimension of the permutation.
  virtual size_type dim() const = 0;

  /// Returns true if \c this is the identity permutation \a I.
  virtual bool is_identity() const = 0;

  /// Prints debug type of information
  virtual std::ostream& output(std::ostream& out) const = 0;

  //@}

  /** @name Vector permutations */
  //@{

  /** \brief Permute a vector <tt>op(P)*x -> y</tt>
   *
   * @param  P_trans  [in] <tt>op(P) = P</tt> for <tt>P_trans == BLAS_Cpp::no_trans</tt> or
   *                  <tt>op(P) = P'</tt> for <tt>P_trans == BLAS_Cpp::trans</tt>.
   * @param  x        [in] Vector.
   * @param  y        [out] Vector.
   *
   * Preconditions:<ul>
   * <li> <tt>y != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>x.space().is_compatible(this->space()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>y->space().is_compatible(this->space()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   */
  virtual void permute( 
    BLAS_Cpp::Transp    P_trans
    ,const Vector       &x
    ,VectorMutable      *y
    ) const = 0;

  /** \brief Permute a vector <tt>op(P)*y -> y</tt>
   *
   * @param  P_trans  [in] <tt>op(P) = P</tt> for <tt>P_trans == BLAS_Cpp::no_trans</tt> or
   *                  <tt>op(P) = P'</tt> for <tt>P_trans == BLAS_Cpp::trans</tt>.
   * @param  y        [in/out] Vector.
   *
   * Preconditions:<ul>
   * <li> <tt>y != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>y->space().is_compatible(this->space()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   */
  virtual void permute( 
    BLAS_Cpp::Transp    P_trans
    ,VectorMutable      *y
    ) const = 0;

  //@}

}; // end class Permutation

} // end namespace AbstractLinAlgPack

#endif // ALAP_PERMUTATION_H
