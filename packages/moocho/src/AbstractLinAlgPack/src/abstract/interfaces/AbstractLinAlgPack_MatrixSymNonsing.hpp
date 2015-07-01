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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_SYM_NONSINGULAR_H
#define ABSTRACT_LINALG_PACK_MATRIX_SYM_NONSINGULAR_H

#include "AbstractLinAlgPack_MatrixNonsing.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all polymorphic symmetrix nonsingular matrices that
 * can be used to solve for linear systems relatively efficently.
 *
 * This interface defines a single addition method to those found in \c MatrixNonsing:
 *
 * <tt>symwo_lhs = alpha * op(mwo) * inv(M) * op(mwo)'</tt><br>
 *
 * The reason that this method could not be defined in the \c MatrixNonsing interface
 * is that the lhs matrix matrix argument \c symwo_lhs is only guaranteed to be
 * symmetric if the rhs matrix argument \c M (which is \c this matrix) is guaranteed
 * to be symmetric.  Since a \c MatrixNonsing matrix object may be unsymmetric, it
 * can not implement this operation, only a symmetric nonsingular matrix can.
 *
 * Any symmetric nonsingular matrix abstraction that can be used to solve for nonlinear
 * systems should also be able to support the \c MatrixSymOp interface.
 * Therefore, this interface is more of an implementation artifact than
 * an a legitimate domain abstraction.  However, some symmetric linear solvers that
 * can implement this interface, can not easily implement the <tt>%MatrixSymOp</tt>
 * interface and therefore this interface is justified.  A general client should never
 * use this interface directly.  Instead, the combined interface \c MatrixSymOpNonsing
 * should be used with fully formed symmetric matrix abstractions.
 *
 * Clients should use the \ref MatrixSymNonsingular_func_grp "provided non-member functions"
 * to call the methods and not the methods themselves.
 */
class MatrixSymNonsing
  : public virtual MatrixNonsing
{
public:

  /** @name Public types */
  //@{

#ifndef DOXYGEN_COMPILE
  /** \brief . */
  typedef Teuchos::RCP<const MatrixSymNonsing>    mat_msns_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<MatrixSymNonsing>          mat_msns_mut_ptr_t;
#endif
  /** \brief . */
  enum EMatrixDummyArg { DUMMY_ARG };

  //@}

  /** @name Friends */
  //@{

  /** \brief . */
  friend
  void M_StMtInvMtM(
    MatrixSymOp* sym_gms_lhs, value_type alpha
    ,const MatrixOp& mwo
    ,BLAS_Cpp::Transp mwo_trans, const MatrixSymNonsing& mswof
    ,EMatrixDummyArg mwo_rhs
     );

  //@}

  /** @name Clone */
  //@{

  /** \brief Clone the non-const matrix object (if supported).
   *
   * The default implementation returns NULL which is perfectly acceptable.
   * A matrix object is not required to return a non-NULL value but almost
   * every good matrix implementation will.
   */
  virtual mat_msns_mut_ptr_t clone_msns();

  /** \brief Clone the const matrix object (if supported).
   *
   * The behavior of this method is the same as for the non-const version
   * above except it returns a smart pointer to a const matrix object.
   */
  virtual mat_msns_ptr_t clone_msns() const;

  //@}

protected:

  /** @name Level-3 */
  //@{

  /** \brief symwo_lhs = alpha * op(mwo) * inv(M) * op(mwo)'.
   *
   * The default implementation is based on the operation M_StInvMtM(...)
   * assuming that this #M# is a symmetric matrix.  For an efficient implementation
   * (for this = L*L' for instance) the subclass may want to override this function.
   */
  virtual void M_StMtInvMtM(
    MatrixSymOp* symwo_lhs, value_type alpha
    ,const MatrixOp& mwo, BLAS_Cpp::Transp mwo_trans
    ,EMatrixDummyArg
    ) const;

  //@}

public:

  /** Overridden from MatrixNonsing */
  //@{
  /// Returns <tt>this->clone_msns()</tt>.
  mat_mns_mut_ptr_t clone_mns();
  /// Returns <tt>this->clone_msns()</tt>.
  mat_mns_ptr_t clone_mns() const;
  //@}

};	// end class MatrixSymNonsing

// ////////////////////////////////////////////////////////////////////////////////////////////////
/** \defgroup MatrixSymNonsingular_func_grp MatrixSymNonsing non-member functions that call virtual functions.
  *
  * These allow nonmember functions to act like virtual functions.
  */
//@{

inline
/// sym_gms_lhs = alpha * op(mwo) * inv(mswof) * op(mwo)'
void M_StMtInvMtM(
  MatrixSymOp* sym_gms_lhs, value_type alpha
  ,const MatrixOp& mwo
  ,BLAS_Cpp::Transp mwo_trans, const MatrixSymNonsing& mswof
  ,MatrixSymNonsing::EMatrixDummyArg mwo_rhs
  )
{
  mswof.M_StMtInvMtM(sym_gms_lhs,alpha,mwo,mwo_trans,mwo_rhs);
}

//@}

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_SYM_NONSINGULAR_H
