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

#ifndef COO_MATRIX_PARTITION_VIEW_SUBCLASS_H
#define COO_MATRIX_PARTITION_VIEW_SUBCLASS_H

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_COOMatrixWithPartitionedView.hpp"

namespace AbstractLinAlgPack {

// Could not derive this class form MatrixWithOpConcreteEncap because the assignment
// operator is not defined for the partition class.

/** \brief Implementation of MatrixOp abstract interface for
  * COOMatrixWithPartitionedView::partition_type.
  */
class COOMatrixPartitionViewSubclass : public MatrixOp
{
public:

  /** \brief . */
  typedef COOMatrixWithPartitionedView::partition_type	M;

  /** \brief . */
  COOMatrixPartitionViewSubclass()
    : trans_(BLAS_Cpp::no_trans)
  {}

  /** \brief . */
  COOMatrixPartitionViewSubclass(BLAS_Cpp::Transp trans)
    : trans_(trans)
  {}

  /** \brief . */
  COOMatrixPartitionViewSubclass(const M& m)
    : m_(m), trans_(BLAS_Cpp::no_trans)
  {}

  /** \brief . */
  COOMatrixPartitionViewSubclass(const M& m, BLAS_Cpp::Transp trans)
    : m_(m), trans_(trans)
  {}

  /** \brief . */
  void set_trans(BLAS_Cpp::Transp trans) {
    trans_ = trans;
  }

  // /////////////////////////////////////////////////////
  /** @name Representation access */
  //@{

  /// Get the underlying M object
  M& m() {
    return m_;
  }

  /** \brief . */
  const M& m() const {
    return m_;
  }

  //		end Representation access
  //@}

  // /////////////////////////////////////////////////////
  // Overridden from Matrix

  /** \brief . */
  size_type rows() const;

  /** \brief . */
  size_type cols() const;

  // /////////////////////////////////////////////////////
  // Overridden from MatrixOp

  /** \brief . */
  MatrixOp& operator=(const MatrixOp& m);

  // /////////////////////////////////////////////////////
  /** @name Level-1 BLAS */
  //@{

  /// (1) gms_lhs += alpha * op(M_rhs) (BLAS xAXPY)
  void Mp_StM(DMatrixSlice* gms_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const;

  //		end Level-1 BLAS
  //@}

  // ////////////////////////////////////////////////////
  /** @name Level-2 BLAS */
  //@{

  /// (2) vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;

  /// (3) vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;

  /// (4) result = vs_rhs1' * op(M_rhs2) * vs_rhs3
  value_type transVtMtV(const DVectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const DVectorSlice& vs_rhs3) const;

  /// (5) result = sv_rhs1' * op(M_rhs2) * sv_rhs3
  value_type transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
    , const SpVectorSlice& sv_rhs3) const;

  //		end Level-2 BLAS
  //@}

  // ////////////////////////////////////////////////////
  /** @name Level-3 BLAS */
  //@{

  /// (6) gms_lhs = alpha * op(M_rhs1) * op(gms_rhs2) + beta * gms_lhs (right) (xGEMM)
  void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha
    , BLAS_Cpp::Transp trans_rhs1, const DMatrixSlice& gms_rhs2
    , BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

  /// (7) gms_lhs = alpha * op(gms_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xGEMM)
  void Mp_StMtM(DMatrixSlice* gms_lhs, value_type alpha, const DMatrixSlice& gms_rhs1
    , BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

  //		end Level-3 BLAS
  //@}

private:
  M m_;
  BLAS_Cpp::Transp trans_;	// for how the matrix if viewed as.

  BLAS_Cpp::Transp op(BLAS_Cpp::Transp trans) const {
    using BLAS_Cpp::trans_not;
    return trans_ == BLAS_Cpp::no_trans ? trans : trans_not(trans);
  }

};	// end class COOMatrixPartitionViewSubclass

}	// end namespace AbstractLinAlgPack 

#endif	// COO_MATRIX_PARTITION_VIEW_SUBCLASS_H
