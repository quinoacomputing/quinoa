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

#ifndef MATRIX_WITH_OP_NONSINGULAR_AGGR_H
#define MATRIX_WITH_OP_NONSINGULAR_AGGR_H

#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"

namespace AbstractLinAlgPack {

/** \brief Aggregate matrix class pulling together a \c MatrixOp object and a
 * \c MatrixNonsing object into a unified matrix object.
 *
 * ToDo: Finish documentation!
 */
class MatrixOpNonsingAggr
  : virtual public MatrixOpNonsing
{
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<const MatrixOp>        mwo_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<const MatrixNonsing>   mns_ptr_t;

  //@}

  /** @name Constructors / initializers */
  //@{

  /// Construct to uninitialized
  MatrixOpNonsingAggr();

  /** \brief Calls <tt>this->initialize()</tt>.
   */
  MatrixOpNonsingAggr(
    const mwo_ptr_t       &mwo
    ,BLAS_Cpp::Transp     mwo_trans
    ,const mns_ptr_t      &mns
    ,BLAS_Cpp::Transp     mns_trans
    );

  /** \brief Initialize.
   *
   * ToDo: Finish documentation.
   */
  void initialize(
    const mwo_ptr_t       &mwo
    ,BLAS_Cpp::Transp     mwo_trans
    ,const mns_ptr_t      &mns
    ,BLAS_Cpp::Transp     mns_trans
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
  const mwo_ptr_t& mwo() const;
  /** \brief . */
  BLAS_Cpp::Transp mwo_trans() const;
  /** \brief . */
  const mns_ptr_t& mns() const;
  /** \brief . */
  BLAS_Cpp::Transp mns_trans() const;

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

  /** @name Overridden from MatrixNonsing */
  //@{

  /** \brief . */
  void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2) const;
  /** \brief . */
  void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2) const;
  /** \brief . */
  value_type transVtInvMtV(
    const Vector& v_rhs1
    ,BLAS_Cpp::Transp trans_rhs2, const Vector& v_rhs3) const;
  /** \brief . */
  value_type transVtInvMtV(
    const SpVectorSlice& sv_rhs1
    ,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const;
  /** \brief . */
  void M_StInvMtM(
    MatrixOp* m_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
    ) const;
  /** \brief . */
  void M_StMtInvM(
    MatrixOp* m_lhs, value_type alpha
    ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2
    ) const;

  //@}

private:

#ifdef DOXYGEN_COMPILE
  MatrixOp         *mwo;
  MatrixNonsing    *mns;
#else
  mwo_ptr_t            mwo_;
  BLAS_Cpp::Transp     mwo_trans_;
  mns_ptr_t            mns_;
  BLAS_Cpp::Transp     mns_trans_;
#endif

}; // end class MatrixOpNonsingAggr

// ////////////////////////////////////
// Inline members

inline
const MatrixOpNonsingAggr::mwo_ptr_t&
MatrixOpNonsingAggr::mwo() const
{
  return mwo_;
}

inline
BLAS_Cpp::Transp MatrixOpNonsingAggr::mwo_trans() const
{
  return mwo_trans_;
}

inline
const MatrixOpNonsingAggr::mns_ptr_t&
MatrixOpNonsingAggr::mns() const
{
  return mns_;
}

inline
BLAS_Cpp::Transp MatrixOpNonsingAggr::mns_trans() const
{
  return mns_trans_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_WITH_OP_NONSINGULAR_AGGR_H
