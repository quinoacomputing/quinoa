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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_SUB_VIEW_H
#define ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_SUB_VIEW_H

#include <iosfwd>

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "Teuchos_RCP.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"

namespace AbstractLinAlgPack {

/** \brief Standard subclass for representing a sub, possibly transposed, view of a matrix
 * 
 * The matrix \c M_view represented by \c this is:
 \verbatim
 
 M_view = op(M_full(rng_rows,rng_cols))
 \endverbatim
 *
 * ToDo: Finish Documentation!
 */
class MatrixOpSubView : public virtual MatrixOp {
public:

  /** \brief . */
  using MatrixOp::syrk;
  
  /** \brief . */
  typedef Teuchos::RCP<MatrixOp>   mat_ptr_t;

  /** @name Constructors/initalizers */
  //@{

  /** \brief Calls <tt>this->initialize(...)</tt>
   */
  MatrixOpSubView(
    const mat_ptr_t&   M_full    = Teuchos::null
    ,const Range1D&    rng_rows  = Range1D()
    ,const Range1D&    rng_cols  = Range1D()
    ,BLAS_Cpp::Transp   M_trans  = BLAS_Cpp::no_trans
    );
    
  /** \brief Initialize the view of a matrix.
   *
   * @param  M_full   [in] Smart pointer for the matrix to provide a view of.
   *                  It is allowed for <tt>M_full.get() == NULL</tt> in which case
   *                  \c this will become uninitialized and none of the rest of the
   *                  arguments matter and any value will do (i.e. the default values).
   * @param  rng_rows [in] Range in the rows of <tt>*M_full</tt> that \c this will represent.
   *                  Only significant if <tt>M_full.get() != NULL</tt>.
   * @param  rng_cols [in] Range in the columns of <tt>*M_full</tt> that \c this will represent.
   *                  Only significant if <tt>M_full.get() != NULL</tt>.
   * @param  M_trans  [in] If <tt>M_trans == no_trans</tt> then \c this will represent
   *                  <tt>M_full(rng_rows,rng_cols)</tt>.  If If <tt>M_trans == trans</tt>
   *                  then \c this will represent the transpose <tt>M_full(rng_rows,rng_cols)'</tt>.
   *
   * Preconditions:<ul>
   * <li>[<tt>M_full.get()!=NULL && !rng_rows.full_range()</tt>]
   *     <tt>rng_rows.ubound() <= M_full->rows()</tt> (throw <tt>std::invalid_argument</tt>)
   * <li>[<tt>M_full.get()!=NULL && !rng_cols.full_range()</tt>]
   *     <tt>rng_cols.ubound() <= M_full->cols()</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->M_full_ptr().get() == M_full.get()</tt>
   * <li>[<tt>M_full.get()!=NULL</tt>] <tt>&this->M_full() == M_full.get()</tt>
   * <li>[<tt>M_full.get()!=NULL && (rng_rows.full_range() || (rng_rows.lbound() == 1 && rng_rows.ubound() == M_full->rows()))</tt>]
   *     <tt>this->rng_rows().full_range() == true</tt>
   * <li>[<tt>M_full.get()!=NULL && (rng_cols.full_range() || (rng_cols.lbound() == 1 && rng_cols.ubound() == M_full->cols()))</tt>]
   *     <tt>this->rng_cols().full_range() == true</tt>
   * <li>[<tt>M_full.get()!=NULL</tt>] <tt>this->M_trans == M_trans</tt>
   * <li>[<tt>M_full.get()==NULL</tt>] <tt>this->rng_rows() == Range1D::Invalid</tt>
   * <li>[<tt>M_full.get()==NULL</tt>] <tt>this->rng_cols() == Range1D::Invalid</tt>
   * </ul>
   */
  void initialize(
    const mat_ptr_t&   M_full
    ,const Range1D&    rng_rows  = Range1D()
    ,const Range1D&    rng_cols  = Range1D()
    ,BLAS_Cpp::Transp  M_trans = BLAS_Cpp::no_trans
    );

  //@}

  /** @name Representation access */
  //@{

  /** \brief . */
  const mat_ptr_t& M_full_ptr();
  /** \brief . */
  MatrixOp& M_full();
  /** \brief . */
  const MatrixOp& M_full() const;
  /** \brief . */
  Range1D rng_rows() const;
  /** \brief . */
  Range1D rng_cols() const;
  /** \brief . */
  BLAS_Cpp::Transp M_trans();

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
  void zero_out();
  /** \brief . */
  void Mt_S( value_type alpha );
  /** \brief . */
  MatrixOp& operator=(const MatrixOp& M);
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
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
  bool Mp_StM(
    value_type alpha,const MatrixOp& M_rhs, BLAS_Cpp::Transp trans_rhs);
  /** \brief . */
  bool Mp_StMtP(
    value_type alpha
    ,const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
    ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    );
  /** \brief . */
  bool Mp_StPtM(
    value_type alpha
    ,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    ,const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
    );
  /** \brief . */
  bool Mp_StPtMtP(
    value_type alpha
    ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    ,const MatrixOp& M_rhs, BLAS_Cpp::Transp M_trans
    ,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
    );
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
  bool Mp_StMtM(
    value_type alpha
    ,const MatrixOp& mvw_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
    ,value_type beta );
  /** \brief . */
  bool syrk(
     BLAS_Cpp::Transp M_trans, value_type alpha
    , value_type beta, MatrixSymOp* sym_lhs ) const;
  
  //@}

private:
  
#ifdef DOXYGEN_COMPILE	
  MatrixOp              *M_full;
  Range1D                   rng_rows;
  Range1D                   rng_cols;
#else
  mat_ptr_t                 M_full_;
  Range1D                   rng_rows_;
  Range1D                   rng_cols_;
  BLAS_Cpp::Transp          M_trans_;
  VectorSpace::space_ptr_t  space_cols_;
  VectorSpace::space_ptr_t  space_rows_;
#endif

  //
  void assert_initialized() const;
  
};	// end class MatrixOpSubView

// //////////////////////////////////
// Inline members

inline
const MatrixOpSubView::mat_ptr_t&
MatrixOpSubView::M_full_ptr()
{
  return M_full_;
}

inline
MatrixOp& MatrixOpSubView::M_full()
{
  return *M_full_;
}

inline
const MatrixOp& MatrixOpSubView::M_full() const
{
  return *M_full_;
}

inline
Range1D MatrixOpSubView::rng_rows() const
{
  return rng_rows_;
}

inline
Range1D MatrixOpSubView::rng_cols() const
{
  return rng_rows_;
}

inline
BLAS_Cpp::Transp MatrixOpSubView::M_trans()
{
  return M_trans_;
}

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_SUB_VIEW_H
