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
//
#ifndef ABSTRACT_LINALG_PACK_MATRIX_NONSINGULAR_H
#define ABSTRACT_LINALG_PACK_MATRIX_NONSINGULAR_H

#include "AbstractLinAlgPack_MatrixBase.hpp"
#include "Teuchos_RCP.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all nonsingular polymorphic matrices that can solve
 * for linear system with but it may not be convienent to compute matrix vector
 * products {abstract}.
 * 
 * The operations supported are:
 *
 * Level-2 BLAS
 *
 * <tt>v_lhs	= inv(op(M_rhs1)) * vs_rhs2</tt><br>
 * <tt>v_lhs	= inv(op(M_rhs1)) * sv_rhs2</tt><br>
 * <tt>result	= v_rhs1' * inv(op(M_rhs2)) * v_rhs3</tt><br>
 * <tt>result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3</tt><br>
 * 
 * Level-3 BLAS
 *
 * <tt>m_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)</tt><br>
 * <tt>m_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)</tt><br>
 *
 * For the solve operations, the lhs and rhs arguments may not be the same
 * in general so don't assume that you can alias the lhs with the rhs and
 * get correct results.
 *
 * Any nonsingular matrix abstraction that can be used to solve for nonlinear
 * systems should also be able to support the \c MatrixOp interface.
 * Therefore, this interface is more of an implementation artifact than
 * an a legitimate domain abstraction.  However, some linear solvers that
 * can implement this interface, can not easily implement the <tt>%MatrixOp</tt>
 * interface and therefore this interface is justified.  A general client should never
 * use this interface directly.  Instead, the combined interface \c MatrixOpNonsing
 * should be used with fully formed matrix abstractions.
 *
 * All these Level-2 and Level-3 BLAS operations have default implementations
 * based on the Level-2 BLAS operations:
 *
 * <tt>v_lhs = inv(op(M_rhs1)) * vs_rhs2</tt><br>
 *
 * which allows for fast prototyping of new matrix subclasses.
 *
 * The member functions should not be called directly but instead through
 * the \ref MatrixNonsingular_funcs_grp "provided non-member functions".
 *
 * The multiple dispatch approach taken in <tt>MatrixOp</tt> is not taken
 * in this interface.  This is because it is considered here that the
 * nonsingular matrix takes procedence of a general matrix arguemnt and
 * we can not expect a general matrix to know how to solve for a linear
 * system with some other nonsigular matrix.
 */
class MatrixNonsing : public virtual MatrixBase {
public:

  /** @name Friends */
  //@{

  /** \brief . */
  friend
  void V_InvMtV(
    VectorMutable* v_lhs, const MatrixNonsing& M_rhs1
    ,BLAS_Cpp::Transp trans_rhs1, const Vector& v_rhs2);
  /** \brief . */
  friend
  void V_InvMtV(
    VectorMutable* v_lhs, const MatrixNonsing& M_rhs1
    ,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2);
  /** \brief . */
  friend
  value_type transVtInvMtV(
    const Vector& v_rhs1, const MatrixNonsing& M_rhs2
    ,BLAS_Cpp::Transp trans_rhs2, const Vector& v_rhs3);
  /** \brief . */
  friend
  value_type transVtInvMtV(
    const SpVectorSlice& sv_rhs1, const MatrixNonsing& M_rhs2
    ,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3);
  /** \brief . */
  friend
  void M_StInvMtM(
    MatrixOp* m_lhs, value_type alpha
    ,const MatrixNonsing&  M_rhs1,     BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp&       mwo_rhs2,   BLAS_Cpp::Transp trans_rhs2 );
  /** \brief . */
  friend
  void M_StMtInvM(
    MatrixOp* m_lhs, value_type alpha
    ,const MatrixOp&      mwo_rhs1,  BLAS_Cpp::Transp trans_rhs1
    ,const MatrixNonsing& M_rhs2,    BLAS_Cpp::Transp trans_rhs2 );

  //@}

  /** @name Public types */
  //@{

#ifndef DOXYGEN_COMPILE
  /** \brief . */
  typedef Teuchos::RCP<const MatrixNonsing>    mat_mns_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<MatrixNonsing>          mat_mns_mut_ptr_t;
#endif

  /** \brief This exception will be thrown if it turns out at runtime that
   * the matrix is numerically singular.
   */
  class SingularMatrix : public std::logic_error
  {public: SingularMatrix(const std::string& what_arg) : std::logic_error(what_arg) {}};

  //@}

  /** @name Clone */
  //@{

  /** \brief Clone the non-const matrix object (if supported).
   *
   * The default implementation returns NULL which is perfectly acceptable.
   * A matrix object is not required to return a non-NULL value but almost
   * every good matrix implementation will.
   */
  virtual mat_mns_mut_ptr_t clone_mns();

  /** \brief Clone the const matrix object (if supported).
   *
   * The behavior of this method is the same as for the non-const version
   * above except it returns a smart pointer to a const matrix object.
   *
   * The default implementation of this method will call the non-const version
   * and then cast to constant.
   */
  virtual mat_mns_ptr_t clone_mns() const;

  //@}

  /** @name Level-2 BLAS */
  //@{

  /// v_lhs	= inv(op(M_rhs1)) * vs_rhs2
  virtual void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2) const = 0;
  /// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
  virtual void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2) const;
  /// result	= vs_rhs1' * inv(op(M_rhs2)) * vs_rhs3
  virtual value_type transVtInvMtV(
    const Vector& v_rhs1
    ,BLAS_Cpp::Transp trans_rhs2, const Vector& v_rhs3) const;
  /// result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3
  virtual value_type transVtInvMtV(
    const SpVectorSlice& sv_rhs1
    ,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const;

  //		end Level-2 BLAS
  //@}

  /** @name Level-3 BLAS */
  //@{

  /** \brief m_lhs = alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right).
   *
   * The default implemention performs a <tt>dynamic_cast<MultiVectorMutable>(m_lhs)</tt>.
   * If this \c dynamic_cast<> does not return  \c NULL , then this operation is implemented in terms of
   * <tt>this->V_InvMtV()</tt> one row or column at a time.  If this \c dynamic_cast<> returns
   * false, then this default implementation has no choice but to throw an exception
   * (<tt>std::invalid_argument</tt>).
   */
  virtual void M_StInvMtM(
    MatrixOp* m_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
    ) const;
  /** \brief m_lhs = alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left).
   *
   * The default implemention performs a <tt>dynamic_cast<MultiVectorMutable>(m_lhs)</tt>.
   * If this \c dynamic_cast<> does not return  \c NULL , then this operation is implemented in terms of
   * <tt>this->V_InvMtV()</tt> one row or column at a time.  If this \c dynamic_cast<> returns
   * false, then this default implementation has no choice but to throw an exception
   * (<tt>std::invalid_argument</tt>).
   */
  virtual void M_StMtInvM(
    MatrixOp* m_lhs, value_type alpha
    ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2
    ) const;

  //		end Level-3 BLAS
  //@}

};	// end class MatrixNonsing

/** \defgroup MatrixNonsingular_funcs_grp MatrixNonsing inline non-member operation functions
  * that call virtual functions.
  *
  * These allow nonmember functions to act like virtual functions
  * and thereby allow the same syntax as in DenseLinAlgPack.
  */
//@{

/** @name Level-2 BLAS */
//@ {

/// v_lhs	= inv(op(M_rhs1)) * v_rhs2
inline void V_InvMtV(
  VectorMutable* v_lhs, const MatrixNonsing& M_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const Vector& v_rhs2)
{
  M_rhs1.V_InvMtV(v_lhs,trans_rhs1,v_rhs2);
}

/// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
inline void V_InvMtV(
  VectorMutable* v_lhs, const MatrixNonsing& M_rhs1
  ,BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2)
{
  M_rhs1.V_InvMtV(v_lhs,trans_rhs1,sv_rhs2);
}

/// result	= v_rhs1' * inv(op(M_rhs2)) * v_rhs3
inline value_type transVtInvMtV(
  const Vector& v_rhs1, const MatrixNonsing& M_rhs2
  ,BLAS_Cpp::Transp trans_rhs2, const Vector& v_rhs3)
{
  return M_rhs2.transVtInvMtV(v_rhs1,trans_rhs2,v_rhs3);
}

/// result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3
inline value_type transVtInvMtV(
  const SpVectorSlice& sv_rhs1, const MatrixNonsing& M_rhs2
  ,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3)
{
  return M_rhs2.transVtInvMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

//		end Level-2 BLAS
//@ }

/** @name Level-3 BLAS */
//@ {

/// m_lhs	= alpha * inv(op(mwo_rhs1)) * op(mwo_rhs2) (right)
inline void M_StInvMtM(
  MatrixOp* m_lhs, value_type alpha
  ,const MatrixNonsing&  M_rhs1,     BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOp&       mwo_rhs2,   BLAS_Cpp::Transp trans_rhs2 )
{
  M_rhs1.M_StInvMtM(m_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2);
}

/// m_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
  MatrixOp* m_lhs, value_type alpha
  ,const MatrixOp&      mwo_rhs1,  BLAS_Cpp::Transp trans_rhs1
  ,const MatrixNonsing& M_rhs2,    BLAS_Cpp::Transp trans_rhs2 )
{
  M_rhs2.M_StMtInvM(m_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2);
}

//		end Level-3 BLAS
//@ }

//		end Inline non-member operation functions
//@}

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_NONSINGULAR_H
