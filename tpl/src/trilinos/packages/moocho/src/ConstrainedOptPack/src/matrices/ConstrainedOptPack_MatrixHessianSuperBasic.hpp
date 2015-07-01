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

#ifndef MATRIX_HESSIAN_SUPER_BASIC_H
#define MATRIX_HESSIAN_SUPER_BASIC_H

#include <vector>

#include "ConstrainedOptPack/src/ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack/src/MatrixSymWithOpFactorized.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "Miref_count_ptr.h"

namespace ConstrainedOptPack {

/** \brief Matrix class that represents a hessian matrix where only the super
 * submatrix for the super basic variables need be nonsingular.
 *
 * Given a Hessian matrix #B# and a partitioning of the variables
 * #Q = [ Q_R  Q_X ]# into free (superbasic) #Q_R'*x# and fixed (nonbasic)
 * #Q_X'*x# variables, this class represents the following matrix:
 \begin{verbatim}
  [n,n] == size(B)
  [n,n] == size(Q), Q is an othogonal permutation matrix (i.e. Q*Q' = Q'*Q = I)
  [n,n_R] == size(Q_R)
  [n,n_X] == size(Q_X)

  B = Q*Q'*B*Q*Q' = [ Q_R  Q_X ] * [ Q_R'*B*Q_R  Q_R'*B*Q_X ] * [ Q_R' ]
                                   [ Q_X'*B*Q_R  Q_X'*B*Q_X ]   [ Q_X' ]

    = [ Q_R  Q_X ] * [   B_RR     op(B_RX) ] * [ Q_R' ]
                     [ op(B_RX')   B_XX    ]   [ Q_X' ]

    = Q_R*B_RR*Q_R' + Q_R*op(B_RX)*Q_X' + Q_X*op(B_RX')*Q_R + Q_X*B_XX*Q_X'
 \end{verbatim}
 * Above, we allow the prepresentation of #op(B_RX) = Q_R'*B*Q_X# to be
 * transposed to allow for more flexibility.  Since #B_RR# and #B_XX# are
 * symmetric, we do not need to worry about transpose or not.  For this class
 * the matrix #B_RR# is required to be symmetric and nonsingular
 * (#MatrixSymWithOpFactorized# interface), but not necessarily positive definite.
 * This is the condition necessary for the Hessian when projected into the
 * active constraints at the solution for an NLP.
 * The other matrices hold no other special properties other than #B_XX# being
 * symmetric of course.
 *
 * The default compiler generated constructors are allowed and initialize the
 * matrix to uninitialized by default.
 */
class MatrixHessianSuperBasic
  : public virtual MatrixSymOp
{
public:

  /** \brief . */
  typedef Teuchos::RCP<const MatrixSymWithOpFactorized>
    B_RR_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const MatrixOp>
    B_RX_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const MatrixSymOp>
    B_XX_ptr_t;
  /** \brief . */
  typedef std::vector<EBounds>
    bnd_fixed_t;

  /** \brief Constructs to uninitialized.
   */
  MatrixHessianSuperBasic();

  /** \brief Initialize the matrix.
   *
   * Preconditions:\begin{itemize}
   * \item [#i_x_free != NULL#] #0 < i_x_free[l-1] <= n, l = 1...n_R# (throw ???)
   * \item [#i_x_fixed != NULL#]#0 < i_x_fixed[l-1] <= n, l = 1...n-n_R# (throw ???)
   * \item [#i_x_free != NULL && i_x_fixed != NULL#]
   *       #i_x_free[l-1] != i_x_fixed[p-1], l = 1...n_R, p = 1...n-n_X# (throw ???)
   * \item [#n_R > 0#] #B_RR_ptr.get() != NULL && B_RR_ptr->rows() == n_R && B_RR_ptr->cols() == n_R#
   *       (throw #std::invalid_argument#)
   * \item [#n_R == 0#] #B_RR_ptr.get() == NULL# (throw #std::invalid_argument#)
   * \item [#n_R < n && B_RX_ptr.get() != NULL && B_RX_trans == no_trans#]
   *       #B_RX_ptr->rows() == n_R && B_RX_ptr->cols() == n-n_R# (throw #std::invalid_argument#)
   * \item [#n_R < n && B_RX_ptr.get() != NULL && B_RX_trans == trans#]
   *       #B_RX_ptr->rows() == n-n_R && B_RX_ptr->cols() == n_R# (throw #std::invalid_argument#)
   * \item [#n_R == n#] #B_RX_ptr.get() == NULL# (throw ##std::invalid_argument#)
   * \item [#n_R < n#] #B_XX_ptr.get() != NULL && B_XX_ptr->rows() == n-n_R && B_XX_ptr->cols() == n-n_R#
   *       (throw #std::invalid_argument#)
   * \item [#n_R == n#] #B_XX_ptr.get() == NULL# (throw ##std::invalid_argument#)
   * \end{itemize}
   *
   * @param  n    [in] number of variables (used for consistency checking)
   * @param  n_R  [in] number of initially free variables (used for consistency checking)
   * @param  i_x_free
   *              [in] array (size #n_R#): #i_x_free[l-1], l = 1...n_R# defines
   *                   the matrix #Q_R# as:\\
   *                   #Q_R(:,l) = e(i_x_free[l-1]), l = 1...n_R#\\
   *                   The ordering of these indices is significant.
   *                   If #n == n_R# then #i_x_free == NULL# is allowed in which case
   *                   it is assumed to be identity.  If #n_R == 0# then of course
   *                   #i_x_free == NULL# is allowed.
   * @param  i_x_fixed
   *              [in] array (size #n_X = n - n_R#):
   *                   #i_x_fixed[l-1], l = 1...n_X# defines the matrix #Q_X# as:\\
   *                   #Q_X(:,l) = e(i_x_fixed[l-1]), l = 1...n_X#\\
   *                   The ordering of these indices is significant.
   *                   If #n_R==0# then #i_x_fixed == NULL# is allowed in which case
   *                   it is assumed to be identity.  If #n_R == n# then of course
   *                   #i_x_fixed == NULL# is allowed.
   * @param  bnd_fixed
   *              [in] array (size #n_X#):
   *                   #bnd_fixed[l-1], l = 1...n_X# defines what bound the variables
   *                   in #i_x_fixed[l-1], l = 1...n_X# are fixed at: #LOWER#, #UPPER#
   *                   or #EQUALITY#.  If #n_R == n# then of course
   *                   #bnd_fixed == NULL# is allowed.
   * @param  B_RR_ptr
   *              [in] Smart pointer to matrix #B_RR# (size #n_R x n_R#) for the
   *                   free (super basic) variables. #B_RR_ptr.get() != NULL#
   *                   must be true if #n_R > # or an exception will be thrown.
   *                   if #n_R == 0# then #B_RR_ptr.get() == NULL# may be true.
   * @param  B_RX_ptr
   *              [in] Smart pointer to matrix #B_RX# (size #n_R x n_X#
   *                   if #B_RX_trans==no_trans#  or #n_X x n_R# if #B_RX_trans==trans#)
   *                   for the  cross terms of free (super basic) and fixed (nonbasic)
   *                   variables.  It is allowed for #B_RX_ptr.get() == NULL#.
   * @param  B_RX_trans
   *              [in] Determines if op(B_RX) = B_RX (#no_trans#) or op(B_RX) = B_RX'
   *                   (#trans#).  Ignored if #n_R == n#.
   * @param  B_XX_ptr
   *              [in] Smart pointer to matrix B_XX (size #n_X x n_X#) for the
   *                   fixed (nonbasic) variables. #B_XX_ptr.get() != NULL#
   *                   must be true if #n_R < n# or an exception will be thrown.
   *                   If #n_R == n# then #B_XX_ptr.get() == NULL# may be true.
   */
  virtual void initialize(
    size_type            n
    ,size_type           n_R
    ,const size_type     i_x_free[]
    ,const size_type     i_x_fixed[]
    ,const EBounds       bnd_fixed[]
    ,const B_RR_ptr_t&   B_RR_ptr
    ,const B_RX_ptr_t&   B_RX_ptr
    ,BLAS_Cpp::Transp    B_RX_trans
    ,const B_XX_ptr_t&   B_XX_ptr
    );

  /** @name Provide access to constituent members */
  //@{

  /** \brief . */
  const GenPermMatrixSlice& Q_R() const;
  /** \brief . */
  const GenPermMatrixSlice& Q_X() const;
  /** \brief . */
  const bnd_fixed_t& bnd_fixed() const;
  /** \brief . */
  const B_RR_ptr_t& B_RR_ptr() const;
  /** \brief . */
  const B_RX_ptr_t& B_RX_ptr() const;
  /** \brief . */
  BLAS_Cpp::Transp B_RX_trans() const;
  /** \brief . */
  const B_XX_ptr_t& B_XX_ptr() const;

  //@}

  /** @name Overridden from Matrix */
  //@{

  /// 
  size_type rows() const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const DVectorSlice& sv_rhs3, value_type beta) const;
  /** \brief . */
  value_type transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const SpVectorSlice& sv_rhs3) const ;

  //@}

protected:

  /** \brief . */
  void assert_initialized() const;

private:

  // ///////////////////////////////////
  // Private types

  typedef std::vector<size_type>		row_i_t;
  typedef std::vector<size_type>		col_j_t;
  
  // ///////////////////////////////////
  // Private data members

  size_type               n_;
  size_type               n_R_;
  GenPermMatrixSlice		Q_R_;        // Sorted by row
  row_i_t					Q_R_row_i_;
  col_j_t					Q_R_col_j_;
  GenPermMatrixSlice		Q_X_;        // Sorted by row
  row_i_t					Q_X_row_i_;
  col_j_t					Q_X_col_j_;
  bnd_fixed_t             bnd_fixed_;
  B_RR_ptr_t              B_RR_ptr_;
  B_RX_ptr_t              B_RX_ptr_;
  BLAS_Cpp::Transp        B_RX_trans_;
  B_XX_ptr_t              B_XX_ptr_;

}; // end class MatrixHessianSuperBasic

// ////////////////////////////////////////////
// Inline members for MatrixHessianSuperBasic

inline
const GenPermMatrixSlice& MatrixHessianSuperBasic::Q_R() const
{
  assert_initialized();
  return Q_R_;
}

inline
const GenPermMatrixSlice& MatrixHessianSuperBasic::Q_X() const
{
  assert_initialized();
  return Q_X_;
}

inline
const MatrixHessianSuperBasic::bnd_fixed_t&
MatrixHessianSuperBasic::bnd_fixed() const
{
  return bnd_fixed_;
}

inline
const MatrixHessianSuperBasic::B_RR_ptr_t&
MatrixHessianSuperBasic::B_RR_ptr() const
{
  assert_initialized();
  return B_RR_ptr_;
}

inline
const MatrixHessianSuperBasic::B_RX_ptr_t&
MatrixHessianSuperBasic::B_RX_ptr() const{
  assert_initialized();
  return B_RX_ptr_;
}

inline
BLAS_Cpp::Transp MatrixHessianSuperBasic::B_RX_trans() const
{
  assert_initialized();
  return B_RX_trans_;
}

inline
const MatrixHessianSuperBasic::B_XX_ptr_t&
MatrixHessianSuperBasic::B_XX_ptr() const
{
  assert_initialized();
  return B_XX_ptr_;
}

} // end namespace ConstrainedOptPack

#endif // MATRIX_HESSIAN_SUPER_BASIC_H
