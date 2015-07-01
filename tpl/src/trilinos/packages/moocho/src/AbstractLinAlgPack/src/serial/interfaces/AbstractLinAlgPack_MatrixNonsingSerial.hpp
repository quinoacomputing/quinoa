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

#ifndef SLAP_MATRIX_NONSINGULAR_SERIAL_H
#define SLAP_MATRIX_NONSINGULAR_SERIAL_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixNonsing.hpp"

namespace AbstractLinAlgPack {

/** \brief Abstract base class for all <tt>AbstractLinAlgPack::MatrixNonsing</tt> objects
 * implemented in shared memory space.
 *
 * This base class does a mapping from fully abstract linear algebra to shared memory
 * linear algebra.
 *
 * These methods should not be called directly but instead should be called through
 * the line \ref MatrixNonsingularSerial_funcs "non-member functions" that are provided.
 */
class MatrixNonsingSerial
  : virtual public AbstractLinAlgPack::MatrixNonsing // doxygen needs full name
{
public:

  /** @name Level-2 BLAS */
  //@{

  /// v_lhs	= inv(op(M_rhs1)) * vs_rhs2
  virtual void V_InvMtV(
    DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const DVectorSlice& vs_rhs2) const;
  /// vs_lhs	= inv(op(M_rhs1)) * vs_rhs2
  virtual void V_InvMtV(
    DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const DVectorSlice& vs_rhs2) const = 0;
  /// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
  virtual void V_InvMtV(
    DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const SpVectorSlice& sv_rhs2) const;
  /// vs_lhs	= inv(op(M_rhs1)) * sv_rhs2
  virtual void V_InvMtV(
    DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const SpVectorSlice& sv_rhs2) const;
  /// result	= vs_rhs1' * inv(op(M_rhs2)) * vs_rhs3
  virtual value_type transVtInvMtV(
    const DVectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2, const DVectorSlice& vs_rhs3) const;
  /// result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3
  virtual value_type transVtInvMtV(
    const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const;

  //		end Level-2 BLAS
  //@}

  /** @name Level-3 BLAS */
  //@{

  /// gm_lhs	= alpha * inv(op(M_rhs1)) * op(gms_rhs2) (right)
  virtual void M_StInvMtM(
    DMatrix* gm_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const DMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2 ) const;
  /// gms_lhs	= alpha * inv(op(M_rhs1)) * op(gms_rhs2) (right)
  virtual void M_StInvMtM(
    DMatrixSlice* gms_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const DMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2 ) const;
  /// gm_lhs	= alpha * op(gms_rhs1) * inv(op(M_rhs2)) (left)
  virtual void M_StMtInvM(
    DMatrix* gm_lhs, value_type alpha
    ,const DMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2 ) const;
  /// gms_lhs	= alpha * op(gms_rhs1) * inv(op(M_rhs2)) (left)
  virtual void M_StMtInvM(
    DMatrixSlice* gms_lhs, value_type alpha
    ,const DMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2 ) const;
  /// gm_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)
  virtual void M_StInvMtM(
    DMatrix* gm_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOpSerial& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2 ) const;
  /// gms_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)
  virtual void M_StInvMtM(
    DMatrixSlice* gms_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOpSerial& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2 ) const;
  /// gm_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
  virtual void M_StMtInvM(
    DMatrix* gm_lhs, value_type alpha
    ,const MatrixOpSerial& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2 ) const;
  /// gms_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
  virtual void M_StMtInvM(
    DMatrixSlice* gms_lhs, value_type alpha
    ,const MatrixOpSerial& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2 ) const;

  //		end Level-3 BLAS
  //@}

  /** Overridden from MatrixNonsing */
  //@{

  /// v_lhs	= inv(op(M_rhs1)) * vs_rhs2
  void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2) const;
  /// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
  void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const SpVectorSlice& sv_rhs2) const;
  /// result	= vs_rhs1' * inv(op(M_rhs2)) * vs_rhs3
  value_type transVtInvMtV(
    const Vector& v_rhs1
    ,BLAS_Cpp::Transp trans_rhs2
    ,const Vector& v_rhs3) const;
  /// m_lhs = alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right).
  void M_StInvMtM(
    MatrixOp* m_lhs, value_type alpha
    ,BLAS_Cpp::Transp trans_rhs1
    ,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
    ) const;
  /// m_lhs = alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left).
  void M_StMtInvM(
    MatrixOp* m_lhs, value_type alpha
    ,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
    ,BLAS_Cpp::Transp trans_rhs2
    ) const;

  //@}

};	// end class MatrixNonsingSerial

/** \defgroup MatrixNonsingularSerial_funcs MatrixNonsingSerial nonmember inline functions.
 *
 * These nonmember functions allow operations to be called on \c MatrixNonsingSerial objects
 * in similar manner to those in \c DenseLinAlgPack.
 */
//@{

/** @name Level-2 BLAS */
//@{

/// v_lhs	= inv(op(M_rhs1)) * vs_rhs2
inline void V_InvMtV(DVector* v_lhs, const MatrixNonsingSerial& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2)
{
  M_rhs1.V_InvMtV(v_lhs,trans_rhs1,vs_rhs2);
}

/// vs_lhs	= inv(op(M_rhs1)) * vs_rhs2
inline void V_InvMtV(DVectorSlice* vs_lhs, const MatrixNonsingSerial& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2)
{
  M_rhs1.V_InvMtV(vs_lhs,trans_rhs1,vs_rhs2);
}

/// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
inline void V_InvMtV(DVector* v_lhs, const MatrixNonsingSerial& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2)
{
  M_rhs1.V_InvMtV(v_lhs,trans_rhs1,sv_rhs2);
}

/// vs_lhs	= inv(op(M_rhs1)) * sv_rhs2
inline void V_InvMtV(DVectorSlice* vs_lhs, const MatrixNonsingSerial& M_rhs1
  , BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2)
{
  M_rhs1.V_InvMtV(vs_lhs,trans_rhs1,sv_rhs2);
}

/// result	= vs_rhs1' * inv(op(M_rhs2)) * vs_rhs3
inline value_type transVtInvMtV(const DVectorSlice& vs_rhs1, const MatrixNonsingSerial& M_rhs2
  , BLAS_Cpp::Transp trans_rhs2, const DVectorSlice& sv_rhs3)
{
  return M_rhs2.transVtInvMtV(vs_rhs1,trans_rhs2,sv_rhs3);
}

/// result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3
inline value_type transVtInvMtV(const SpVectorSlice& sv_rhs1, const MatrixNonsingSerial& M_rhs2
  , BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3)
{
  return M_rhs2.transVtInvMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

//		end Level-2 BLAS
//@}

/** @name Level-3 BLAS */
//@{

/// gm_lhs	= alpha * inv(op(M_rhs1)) * op(gms_rhs2) (right)
inline void M_StInvMtM(
  DMatrix* gm_lhs, value_type alpha
  ,const MatrixNonsingSerial& M_rhs1,   BLAS_Cpp::Transp trans_rhs1
  ,const DMatrixSlice&          gms_rhs2, BLAS_Cpp::Transp trans_rhs2
  )
{
  M_rhs1.M_StInvMtM(gm_lhs,alpha,trans_rhs1,gms_rhs2,trans_rhs2);
}

/// gms_lhs	= alpha * inv(op(M_rhs1)) * op(gms_rhs2) (right)
inline void M_StInvMtM(
  DMatrixSlice* gms_lhs, value_type alpha
  ,const MatrixNonsingSerial& M_rhs1,   BLAS_Cpp::Transp trans_rhs1
  ,const DMatrixSlice&          gms_rhs2, BLAS_Cpp::Transp trans_rhs2
  )
{
  M_rhs1.M_StInvMtM(gms_lhs,alpha,trans_rhs1,gms_rhs2,trans_rhs2);
}

/// gm_lhs	= alpha * op(gms_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
  DMatrix* gm_lhs, value_type alpha
  ,const DMatrixSlice&          gms_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixNonsingSerial& M_rhs2,   BLAS_Cpp::Transp trans_rhs2
  )
{
  M_rhs2.M_StMtInvM(gm_lhs,alpha,gms_rhs1,trans_rhs1,trans_rhs2);
}

/// gms_lhs	= alpha * op(gms_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
  DMatrixSlice* gms_lhs, value_type alpha
  ,const DMatrixSlice&          gms_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixNonsingSerial& M_rhs2,   BLAS_Cpp::Transp trans_rhs2
  )
{
  M_rhs2.M_StMtInvM(gms_lhs,alpha,gms_rhs1,trans_rhs1,trans_rhs2);
}

/// gm_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)
inline void M_StInvMtM(
  DMatrix* gm_lhs, value_type alpha
  ,const MatrixNonsingSerial& M_rhs1,   BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOpSerial&      mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  )
{
  M_rhs1.M_StInvMtM(gm_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2);
}

/// gms_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)
inline void M_StInvMtM(
  DMatrixSlice* gms_lhs, value_type alpha
  ,const MatrixNonsingSerial& M_rhs1,   BLAS_Cpp::Transp trans_rhs1
  ,const MatrixOpSerial&      mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
  )
{
  M_rhs1.M_StInvMtM(gms_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2);
}

/// gm_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
  DMatrix* gm_lhs, value_type alpha
  ,const MatrixOpSerial&      mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixNonsingSerial& M_rhs2,   BLAS_Cpp::Transp trans_rhs2
  )
{
  M_rhs2.M_StMtInvM(gm_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2);
}

/// gms_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
  DMatrixSlice* gms_lhs, value_type alpha
  ,const MatrixOpSerial&      mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
  ,const MatrixNonsingSerial& M_rhs2,   BLAS_Cpp::Transp trans_rhs2
  )
{
  M_rhs2.M_StMtInvM(gms_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2);
}

//		end Level-3 BLAS
//@}

//		end Inline non-member operation functions
//@}

}	// end namespace AbstractLinAlgPack

#endif	// SLAP_MATRIX_NONSINGULAR_SERIAL_H
