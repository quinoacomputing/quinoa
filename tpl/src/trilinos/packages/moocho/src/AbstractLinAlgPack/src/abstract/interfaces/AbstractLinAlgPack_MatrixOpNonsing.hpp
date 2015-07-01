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

#ifndef ALAP_MATRIX_OP_NONSING_HPP
#define ALAP_MATRIX_OP_NONSING_HPP

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_MatrixNonsing.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all nonsingular polymorphic matrices
 * that can be used to compute matrix-vector products and solve
 * for linear systems efficiently.
 */
class MatrixOpNonsing
  : virtual public MatrixOp
  , virtual public MatrixNonsing
{
public:

  /** @name Public types */
  //@{

#ifndef DOXYGEN_COMPILE
  /** \brief . */
  typedef Teuchos::RCP<const MatrixOpNonsing>    mat_mwons_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<MatrixOpNonsing>          mat_mwons_mut_ptr_t;
#endif

  //@}

  /** @name Clone */
  //@{

  /** \brief Clone the non-const matrix object (if supported).
   *
   * The default implementation returns NULL which is perfectly acceptable.
   * A matrix object is not required to return a non-NULL value but almost
   * every good matrix implementation will.
   */
  virtual mat_mwons_mut_ptr_t clone_mwons();

  /** \brief Clone the const matrix object (if supported).
   *
   * The behavior of this method is the same as for the non-const version
   * above except it returns a smart pointer to a const matrix object.
   */
  virtual mat_mwons_ptr_t clone_mwons() const;

  //@}

  /** @name Condition number estimation */
  //@{

  /** \brief Compute an estimate of the condition number of this matrix.
   *
   * @param  requested_norm_type
   *                    [in] Determines the requested type of norm for the condition number.
   * @param  allow_replacement
   *                    [in] Determines if the requested norm in specified in <tt>norm_type</tt>
   *                    can be replaced with another norm that can be computde by the matrix.
   *
   * @return If a condition number is computed, then <tt>return.value</tt> gives the value of
   * the condition number in the norm of type <tt>return.type</tt>.
   *
   * Postconditions:<ul>
   * <li> If <tt>allow_replacement==true</tt>, the matrix object must return a computed
   *      condition number who's type is given in <tt>return.type</tt>.
   * <li> If <tt>allow_replacement==false</tt> and the underlying matrix object can not compute
   *      condition number for the norm requested in <tt>norm_type</tt>, then a
   *      <tt>MethodNotImplemented</tt> exception will be thrown.  If the matrix object can
   *      an estimate of the condition number for this norm, then <tt>return.type</tt>
   *      will be equal to <tt>requested_norm_type</tt>.
   * </ul>
   *
   * The default implementation of this method uses Algorithm 2.5 in
   * "Applied Numerical Linear Algebra" by James Demmel (1997) to
   * estimate ||inv(M)||1 or ||inv(M)||inf.  The algorithm uses some
   * of the refinements in the referenced algorithm by Highman.
   * This algorithm only requires solves and transposed solves so
   * every nonsingular matrix object can implement this method.  The
   * default arguments for this function will compute an estimate of
   * the condition number and will not thrown an exception.  The
   * default implementation will throw an exception for any other
   * norm type than <tt>requested_norm_type = MAT_NORM_1</tt> or
   * <tt>requested_norm_type = MAT_NORM_INF</tt>.
   */
  const MatNorm calc_cond_num(
    EMatNormType  requested_norm_type = MAT_NORM_1
    ,bool         allow_replacement   = false
    ) const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{
  /// Returns <tt>this->clone_mwons()</tt>.
  mat_mut_ptr_t clone();
  /// Returns <tt>this->clone_mwons()</tt>.
  mat_ptr_t clone() const;
  //@}

  /** @name Overridden from MatrixNonsing */
  //@{
  /// Returns <tt>this->clone_mwons()</tt>.
  mat_mns_mut_ptr_t clone_mns();
  /// Returns <tt>this->clone_mwons()</tt>.
  mat_mns_ptr_t clone_mns() const;
  //@}

  /// Calls operator=(MatrixOp&)
  MatrixOpNonsing& operator=(const MatrixOpNonsing& M)
  { static_cast<MatrixOp*>(this)->operator=(M); return *this; }

}; // end class MatrixOpNonsing

}	// end namespace AbstractLinAlgPack

#endif	// ALAP_MATRIX_OP_NONSING_HPP
