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

#ifndef MATRIX_PERM_AGGR_H
#define MATRIX_PERM_AGGR_H

#include "AbstractLinAlgPack_MatrixOp.hpp"

namespace AbstractLinAlgPack {

/** \brief Aggregate matrix class for a matrix and its permuted view.
 *
 * <tt>mat_perm = row_perm * mat_orig * col_perm'</tt>.
 *
 */
class MatrixPermAggr
  : virtual public MatrixOp
{
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<const Permutation>   perm_ptr_t;

  //@}

  /** @name Constructors / initializers */
  //@{

  /// Construct to uninitialized
  MatrixPermAggr();

  /** \brief Calls <tt>this->initialize()</tt>.
   */
  MatrixPermAggr(
    const mat_ptr_t      &mat_orig
    ,const perm_ptr_t    &row_perm
    ,const perm_ptr_t    &col_perm
    ,const mat_ptr_t     &mat_perm
    );

  /** \brief Initialize.
   *
   * <tt>mat_perm = row_perm' * mat_orig * col_perm</tt>.
   *
   * @param  mat_orig  [in] Smart pointer to original unpermuted matrix.
   * @param  row_perm  [in] Smart pointer to row permutation.  If <tt>row_perm.get() == NULL</tt>
   *                   then the identity permutation is assumed.
   * @param  col_perm  [in] Smart pointer to column permutation.  If <tt>col_perm.get() == NULL</tt>
   *                   then the identity permutation is assumed.
   * @param  mat_perm  [in] Smart pointer to permuted matrix.  It is allowed for
   *                   <tt>mat_perm.get() == NULL</tt> in which case all of the linear algebra
   *                   methods are implemented in terms of \c mat_orig, \c row_perm and \c col_perm.
   *
   * Preconditions:<ul>
   * <li> <tt>mat_perm.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> [<tt>row_perm.get() != NULL</tt>] <tt>mat_orig->space_cols().is_compatible(row_perm->space()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> [<tt>col_perm.get() != NULL</tt>] <tt>mat_orig->space_rows().is_compatible(col_perm->space()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->mat_orig().get() == mat_orig.get()</tt>
   * <li> <tt>this->row_perm().get() == row_perm.get()</tt>
   * <li> <tt>this->col_perm().get() == col_perm.get()</tt>
   * <li> <tt>this->mat_perm().get() == mat_perm.get()</tt>
   * </ul>
   */
  void initialize(
    const mat_ptr_t      &mat_orig
    ,const perm_ptr_t    &row_perm
    ,const perm_ptr_t    &col_perm
    ,const mat_ptr_t     &mat_perm
    );

  /** \brief Set uninitialized.
   *
   * ToDo: Finish documentation.
   */
  void set_uninitialized();

  //@}
  
  /** @name Access */
  //@{

  /** \brief . */
  const mat_ptr_t& mat_orig() const;
  /** \brief . */
  const perm_ptr_t& row_perm() const;
  /** \brief . */
  const perm_ptr_t& col_perm() const;
  /** \brief . */
  const mat_ptr_t& mat_perm() const;

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
  MatrixOp::mat_ptr_t sub_view(const Range1D& row_rng, const Range1D& col_rng) const;
  /** \brief . */
  MatrixOp& operator=(const MatrixOp& M);
  /** \brief . */
  std::ostream& output(std::ostream& out) const;

protected:

  /** \brief . */
  bool Mp_StM(
    MatrixOp* mwo_lhs, value_type alpha
    , BLAS_Cpp::Transp trans_rhs) const;
  /** \brief . */
  bool Mp_StMtP(
    MatrixOp* mwo_lhs, value_type alpha
    , BLAS_Cpp::Transp M_trans
    , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    ) const;
  /** \brief . */
  bool Mp_StPtM(
    MatrixOp* mwo_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    , BLAS_Cpp::Transp M_trans
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
    , const Vector& v_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const Vector& v_rhs3, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const SpVectorSlice& sv_rhs3, value_type beta) const;
  /** \brief . */
  value_type transVtMtV(
    const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const Vector& v_rhs3) const;
  /** \brief . */
  value_type transVtMtV(
    const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const SpVectorSlice& sv_rhs3) const;
  /** \brief . */
  void syr2k(
     BLAS_Cpp::Transp M_trans, value_type alpha
    , const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
    , const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
    , value_type beta, MatrixSymOp* symwo_lhs ) const;
  /** \brief . */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    , BLAS_Cpp::Transp trans_rhs1, const MatrixOp& mwo_rhs2
    , BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;
  /** \brief . */
  bool Mp_StMtM(
    MatrixOp* mwo_lhs, value_type alpha
    , const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    , BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;
  /** \brief . */
  bool syrk(
     BLAS_Cpp::Transp M_trans, value_type alpha
    , value_type beta, MatrixSymOp* sym_lhs ) const;
  
  //@}

private:

#ifdef DOXYGEN_COMPILE
  MatrixOp         *mat_orig;
  Permutation          *row_perm;
  Permutation          *col_perm;
  MatrixOp         *mat_perm;
#else
  mat_ptr_t            mat_orig_;
  perm_ptr_t           row_perm_;
  perm_ptr_t           col_perm_;
  mat_ptr_t            mat_perm_;
#endif

}; // end class MatrixPermAggr

// ////////////////////////////////////
// Inline members

inline
const MatrixOp::mat_ptr_t&
MatrixPermAggr::mat_orig() const
{
  return mat_orig_;
}

inline
const MatrixPermAggr::perm_ptr_t&
MatrixPermAggr::row_perm() const
{
  return row_perm_;
}

inline
const MatrixPermAggr::perm_ptr_t&
MatrixPermAggr::col_perm() const
{
  return col_perm_;
}

inline
const MatrixOp::mat_ptr_t& MatrixPermAggr::mat_perm() const
{
  return mat_perm_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_PERM_AGGR_H
