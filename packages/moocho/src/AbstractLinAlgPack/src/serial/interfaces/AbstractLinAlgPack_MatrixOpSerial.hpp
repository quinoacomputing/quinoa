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

#ifndef SPARSE_LINALG_PACK_MATRIX_WITH_OP_SERIAL_H
#define SPARSE_LINALG_PACK_MATRIX_WITH_OP_SERIAL_H

#include <iosfwd>

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorSpaceSerial.hpp"

namespace AbstractLinAlgPack {

/** \brief Base class for all matrices implemented in a shared memory address space.
 * 
 * This base class inherites from <tt>AbstractLinAlgPack::MatrixOp</tt> and
 * provides a means of transition from purely abstract (i.e. parallel, out-of-core)
 * linear algebra to serial linear algebra (with explicit dense vectors and
 * matrices.
 *
 * This matrix subclass simply uses <tt>dynamic_cast<></tt> to map from the abstract interfaces
 * <tt>AbstractLinAlgPack::Vector</tt>, <tt>AbstractLinAlgPack::VectorMutable</tt>,
 * <tt>AbstractLinAlgPack::MatrixOp</tt> to a concreate types using a few standarized serial
 * interfaces <tt>VectorGetSparse</tt>, <tt>VectorMutableDense</tt> ...
 * If these interfaces are not supported by the vector and matrix arguments, then exceptions may
 * be thrown.
 *
 * Note that every method from \c MatrixOp is not overridden here.  Specifically, the methods
 * <tt>MatrixOp::zero_out()</tt>  and <tt>MatrixOp::sub_view()</tt>
 * are not overridden.  These methods have default implementations that should be sufficient
 * for most uses but subclasses may want to provide specialized implementations anyway.
 *
 * These methods should not be called directly but instead through a set of
 * \ref MatrixWithOpSerial_funcs "provided nonmember functions".
 *
 * ToDo: Finish documentation!
 */
class MatrixOpSerial
  : virtual public AbstractLinAlgPack::MatrixOp // doxygen needs full name
{
public:

  /** \brief . */
  using MatrixOp::Mp_StM;
  /** \brief . */
  using MatrixOp::Mp_StMtP;
  /** \brief . */
  using MatrixOp::Mp_StPtM;
  /** \brief . */
  using MatrixOp::Mp_StPtMtP;
  /** \brief . */
  using MatrixOp::Vp_StMtV;
  /** \brief . */
  using MatrixOp::Mp_StMtM;

  /** @name Level-1 BLAS */
  //@{

  /// gms_lhs += alpha * op(M_rhs) (BLAS xAXPY)
  virtual void Mp_StM(DMatrixSlice* gms_lhs, value_type alpha
    , BLAS_Cpp::Transp trans_rhs) const;

  /// gms_lhs += alpha * op(M_rhs) * op(P_rhs)
  virtual void Mp_StMtP(DMatrixSlice* gms_lhs, value_type alpha
    , BLAS_Cpp::Transp M_trans
    , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    ) const;

  /// gms_lhs += alpha * op(P) * op(M_rhs)
  virtual void Mp_StPtM(DMatrixSlice* gms_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
    , BLAS_Cpp::Transp M_trans
    ) const;

  /// gms_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
  virtual void Mp_StPtMtP(DMatrixSlice* gms_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_trans
    , const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
    ) const;

  //@}

  /** @name Level-2 BLAS */
  //@{

  /// vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
  virtual void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const = 0;

  /// vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
  virtual void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;

  /// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * vs_rhs3 + beta * vs_rhs
  virtual void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const DVectorSlice& vs_rhs3, value_type beta) const;

  /// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * sv_rhs3 + beta * vs_rhs
  virtual void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const SpVectorSlice& sv_rhs3, value_type beta) const;

  /// result = vs_rhs1' * op(M_rhs2) * vs_rhs3
  virtual value_type transVtMtV(const DVectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const DVectorSlice& vs_rhs3) const;

  /** \brief Perform a specialized rank-2k update of a dense symmetric matrix of the form:
   *
   * #sym_lhs += alpha*op(P1')*op(M)*op(P2) + alpha*op(P2')*op(M')*op(P1) + beta*sym_lhs#
   *
   * The reason that this operation is being classified as a level-2 operation is that the
   * total flops should be of O(n^2) and not O(n^2*k).
   *
   * The default implementation is based on #Mp_StMtP(...)# and #Mp_StPtM(...)#.  Of course
   * in situations where this default implemention is inefficient the subclass should
   * override this method.
   */
  virtual void syr2k(
     BLAS_Cpp::Transp M_trans, value_type alpha
    , const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
    , const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
    , value_type beta, DMatrixSliceSym* sym_lhs ) const;

  //@}

  /** @name Level-3 BLAS */
  //@{

  /// gms_lhs = alpha * op(M_rhs1) * op(gms_rhs2) + beta * gms_lhs (right) (xGEMM)
  virtual void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha
    , BLAS_Cpp::Transp trans_rhs1, const DMatrixSlice& gms_rhs2
    , BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

  /// gms_lhs = alpha * op(gms_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xGEMM)
  virtual void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSlice& gms_rhs1
    , BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

  /// gms_lhs = alpha * op(M_rhs1) * op(mwo_rhs2) + beta * gms_lhs (right) (xGEMM)
  virtual void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha
    , BLAS_Cpp::Transp trans_rhs1, const MatrixOpSerial& mwo_rhs2
    , BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

  /// gms_lhs = alpha * op(M_rhs1) * op(sym_rhs2) + beta * gms_lhs (right) (xSYMM)
  virtual void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha
    , BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceSym& sym_rhs2
    , BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

  /// gms_lhs = alpha * op(sym_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xSYMM)
  virtual void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceSym& sym_rhs1
    , BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

  /// gms_lhs = alpha * op(M_rhs1) * op(tri_rhs2) + beta * gms_lhs (right) (xTRMM)
  virtual void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha
    , BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceTri& tri_rhs2
    , BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

  /// gms_lhs = alpha * op(tri_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xTRMM)
  virtual void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceTri& tri_rhs1
    , BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;


  /** \brief Perform a rank-k update of a dense symmetric matrix of the form:
   *
   * #sym_lhs += op(M)*op(M') + beta*sym_lhs#
   *
   * The default implementation is based on #Vp_StMtV(...)#.  Of course
   * in situations where this default implemention is inefficient
   * the subclass should override this method.
   */
  virtual void syrk(
     BLAS_Cpp::Transp M_trans, value_type alpha
    , value_type beta, DMatrixSliceSym* sym_lhs ) const;

  //@}

  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  const VectorSpace& space_rows() const;
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
    VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const Vector& v_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* v_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const Vector& v_rhs3, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* v_lhs, value_type alpha
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
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
    ,value_type beta ) const;
  /** \brief . */
  bool syrk(
    BLAS_Cpp::Transp M_trans, value_type alpha
    ,value_type beta, MatrixSymOp* sym_lhs ) const;

  //@}

private:

  // ////////////////////////////////////
  // Private data members

  mutable VectorSpaceSerial       space_cols_;
  mutable VectorSpaceSerial       space_rows_;

  // ////////////////////////////////////
  // Private member functions

};	// end class MatrixOpSerial

/** \defgroup MatrixWithOpSerial_funcs MatrixOpSerial nonmember inline functions.
 *
 * These nonmember functions allow operations to be called on \c MatrixOpSerial objects
 * in similar manner to those in \c DenseLinAlgPack.
 */
//@{

/** @name Level-1 BLAS */
//@{

/// gms_lhs += alpha * op(M_rhs) (BLAS xAXPY)
inline void Mp_StM(DMatrixSlice* gms_lhs, value_type alpha, const MatrixOpSerial& M_rhs
  , BLAS_Cpp::Transp trans_rhs)
{
  M_rhs.Mp_StM(gms_lhs,alpha,trans_rhs);
}

/// gms_lhs += alpha * op(M_rhs) * op(P_rhs)
inline void Mp_StMtP(DMatrixSlice* gms_lhs, value_type alpha
  , const MatrixOpSerial& M_rhs, BLAS_Cpp::Transp M_trans
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  )
{
  M_rhs.Mp_StMtP(gms_lhs,alpha,M_trans,P_rhs,P_rhs_trans);
}

/// gms_lhs += alpha * op(P) * op(M_rhs)
inline void Mp_StPtM(DMatrixSlice* gms_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
  , const MatrixOpSerial& M_rhs, BLAS_Cpp::Transp M_trans
  )
{
  M_rhs.Mp_StPtM(gms_lhs,alpha,P_rhs,P_rhs_trans,M_trans);
}

/// gms_lhs += alpha * op(P_rhs1) * op(M_rhs) * op(P_rhs2)
inline void Mp_StPtMtP(DMatrixSlice* gms_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , const MatrixOpSerial& M_rhs, BLAS_Cpp::Transp trans_rhs
  , const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
  )
{
  M_rhs.Mp_StPtMtP(gms_lhs,alpha,P_rhs1,P_rhs1_trans,trans_rhs,P_rhs2,P_rhs2_trans);
}

//@}

/** @name Level-2 BLAS */
//@{

/// vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
inline void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, const MatrixOpSerial& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2, value_type beta = 1.0)
{
  M_rhs1.Vp_StMtV(vs_lhs,alpha,trans_rhs1,vs_rhs2,beta);
}

/// vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
inline void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, const MatrixOpSerial& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2, value_type beta = 1.0)
{
  M_rhs1.Vp_StMtV(vs_lhs,alpha,trans_rhs1,sv_rhs2,beta);
}

/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * vs_rhs3 + beta * vs_rhs
inline void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , const MatrixOpSerial& M_rhs2, BLAS_Cpp::Transp M_rhs2_trans
  , const DVectorSlice& vs_rhs3, value_type beta = 1.0) 
{
  M_rhs2.Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,vs_rhs3,beta);
}

/// vs_lhs = alpha * op(P_rhs1) * op(M_rhs2) * sv_rhs3 + beta * vs_rhs
inline void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
  , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
  , const MatrixOpSerial& M_rhs2, BLAS_Cpp::Transp M_rhs2_trans
  , const SpVectorSlice& sv_rhs3, value_type beta = 1.0)
{
  M_rhs2.Vp_StPtMtV(vs_lhs,alpha,P_rhs1,P_rhs1_trans,M_rhs2_trans,sv_rhs3,beta);
}

/// result = vs_rhs1' * op(M_rhs2) * vs_rhs3
inline value_type transVtMtV(const DVectorSlice& vs_rhs1, const MatrixOpSerial& M_rhs2
  , BLAS_Cpp::Transp trans_rhs2, const DVectorSlice& vs_rhs3)
{
  return M_rhs2.transVtMtV(vs_rhs1,trans_rhs2,vs_rhs3);
}

/// result = sv_rhs1' * op(M_rhs2) * sv_rhs3
inline value_type transVtMtV(const SpVectorSlice& sv_rhs1, const MatrixOpSerial& M_rhs2
  , BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3)
{
  return M_rhs2.transVtMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

//@}

/** @name Level-3 BLAS */
//@{

/// gms_lhs = alpha * op(M_rhs1) * op(gms_rhs2) + beta * gms_lhs (right) (xGEMM)
inline void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const MatrixOpSerial& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DMatrixSlice& gms_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta = 1.0)
{
  M_rhs1.Mp_StMtM(gms_lhs,alpha,trans_rhs1,gms_rhs2,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(gms_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xGEMM)
inline void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSlice& gms_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOpSerial& M_rhs2, BLAS_Cpp::Transp trans_rhs2
  , value_type beta = 1.0)
{
  M_rhs2.Mp_StMtM(gms_lhs,alpha,gms_rhs1,trans_rhs1,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(mwo_rhs1) * op(mwo_rhs2) + beta * gms_lhs (right) (xGEMM)
inline void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const MatrixOpSerial& mwo_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOpSerial& gms_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta = 1.0)
{
  mwo_rhs1.Mp_StMtM(gms_lhs,alpha,trans_rhs1,gms_rhs2,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(M_rhs1) * op(sym_rhs2) + beta * gms_lhs (right) (xSYMM)
inline void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const MatrixOpSerial& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceSym& sym_rhs2
  , BLAS_Cpp::Transp trans_rhs2, value_type beta = 1.0)
{
  M_rhs1.Mp_StMtM(gms_lhs,alpha,trans_rhs1,sym_rhs2,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(sym_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xSYMM)
inline void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceSym& sym_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOpSerial& M_rhs2, BLAS_Cpp::Transp trans_rhs2
  , value_type beta = 1.0)
{
  M_rhs2.Mp_StMtM(gms_lhs,alpha,sym_rhs1,trans_rhs1,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(M_rhs1) * op(tri_rhs2) + beta * gms_lhs (right) (xTRMM)
inline void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const MatrixOpSerial& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DMatrixSliceTri& tri_rhs2, BLAS_Cpp::Transp trans_rhs2
  , value_type beta = 1.0)
{
  M_rhs1.Mp_StMtM(gms_lhs,alpha,trans_rhs1,tri_rhs2,trans_rhs2,beta);
}

/// gms_lhs = alpha * op(tri_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xTRMM)
inline void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSliceTri& tri_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const MatrixOpSerial& M_rhs2, BLAS_Cpp::Transp trans_rhs2
  , value_type beta = 1.0)
{
  M_rhs2.Mp_StMtM(gms_lhs,alpha,tri_rhs1,trans_rhs1,trans_rhs2,beta);
}

//@}

//@}

}	// end namespace AbstractLinAlgPack 

#endif	// SPARSE_LINALG_PACK_MATRIX_WITH_OP_SERIAL_H
