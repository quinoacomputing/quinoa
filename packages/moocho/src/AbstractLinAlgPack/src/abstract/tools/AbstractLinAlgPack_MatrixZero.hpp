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

#ifndef ALAP_MATRIX_ZERO_H
#define ALAP_MATRIX_ZERO_H

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Implementation of a matrix with all zeros.
 *
 * This may seem like a silly class but it is helpful in some circumstances.
 * This class needs to be updated whenever methods are added or removed
 * from \c MatrixOp.
 */
class MatrixZero : public MatrixOp {
public:

  /** @name Constructors/initializers */
  //@{

  /// Calls <tt>this->initalize()</tt>
  MatrixZero(
    const VectorSpace::space_ptr_t&     space_cols = Teuchos::null
    ,const VectorSpace::space_ptr_t&    space_rows = Teuchos::null
    );

  /** \brief Initialize (or initialize) given the columns and rows vector spaces.
   *
   * Preconditions:<ul>
   * <li> [<tt>space_cols.get() == NULL</tt>] <tt><tt>space_rows.get() == NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>space_cols.get() != NULL</tt>] <tt><tt>space_rows.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>space_cols.get() == NULL</tt>] <tt>this->rows() == 0</tt>
   * <li> [<tt>space_cols.get() == NULL</tt>] <tt>this->cols() == 0</tt>
   * <li> [<tt>space_cols.get() != NULL</tt>] <tt>&this->space_cols() == space_cols.get()</tt>
   * <li> [<tt>space_cols.get() != NULL</tt>] <tt>&this->space_rows() == space_rows.get()</tt>
   * </ul>
   */
  void initialize(
    const VectorSpace::space_ptr_t&    space_cols
    ,const VectorSpace::space_ptr_t&   space_rows
    );

  //@}

  /** @name Overridden from MatrixBase */
  //@{
  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  size_type nz() const;
  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  const VectorSpace& space_rows() const;
  /** \brief . */
  void zero_out();
  /** \brief . */
  void Mt_S(value_type alpha);
  /** \brief . */
  MatrixOp& operator=(const MatrixOp& M);
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  bool Mp_StM(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs
    ) const;
  /** \brief . */
  bool Mp_StMtP(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp M_trans
    ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    ) const;
  /** \brief . */
  bool Mp_StPtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    ,BLAS_Cpp::Transp M_trans
    ) const;
  /** \brief . */
  bool Mp_StPtMtP(
    MatrixOp* mwo_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    ,BLAS_Cpp::Transp M_trans
    ,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
    ) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* vs_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    ,BLAS_Cpp::Transp M_rhs2_trans
    ,const Vector& v_rhs3, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* vs_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    ,BLAS_Cpp::Transp M_rhs2_trans
    ,const SpVectorSlice& sv_rhs3, value_type beta) const;
  /** \brief . */
  value_type transVtMtV(
    const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
    ,const Vector& v_rhs3) const;
  /** \brief . */
  value_type transVtMtV(
    const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const SpVectorSlice& sv_rhs3) const;
  /** \brief . */
  void syr2k(
     BLAS_Cpp::Transp M_trans, value_type alpha
    ,const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
    ,const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
    ,value_type beta, MatrixSymOp* symwo_lhs ) const;
  /** \brief . */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
    ,value_type beta ) const;
  /** \brief . */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2
    ,value_type beta ) const;
  /** \brief . */
  bool syrk(
     BLAS_Cpp::Transp M_trans, value_type alpha
    ,value_type beta, MatrixSymOp* sym_lhs ) const;

  //@}

private:

  VectorSpace::space_ptr_t  space_cols_;
  VectorSpace::space_ptr_t  space_rows_;

  //
  void assert_initialized() const;

}; // end class MatrixZero

} // end namespace AbstractLinAlgPack

#endif // ALAP_MATRIX_ZERO_H
