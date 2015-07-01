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

#ifndef GEN_MATRIX_AS_TRI_SYM_H
#define GEN_MATRIX_AS_TRI_SYM_H

#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace DenseLinAlgPack {

/* * @name	Packaging arguments for GenMatrixSlices treated as triangular
  *			and symmetric matrices in BLAS-like linear algebra operations.
  *
  */

// @{

// /////////////////////////////////////////////////////////////////////////////////////
/** \brief . */
/* * Aggregate information for a triangular matrix (element-wise) stored in a DMatrix.
  *
  * This is the type to be used as lhs and rhs arguments in element-wise
  * linear algebra operations like assignment and binary arithmetic.
  */
class DMatrixSliceTriEle {
public:
  /** \brief . */
  DMatrixSliceTriEle(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
    : gms_(const_cast<DMatrixSlice&>(gms)), uplo_(uplo)
  {
    #ifdef LINALGPACK_CHECK_RHS_SIZES
      assert_gms_square(gms);
    #endif
  }
  /** \brief . */
  size_type rows() const {
    return gms_.rows();
  }
  /** \brief . */
  size_type cols() const {
    return gms_.cols();
  }
  /** \brief . */
  DMatrixSlice& gms() {
    return gms_;
  }
  /** \brief . */
  const DMatrixSlice& gms() const {
    return gms_;
  }
  /** \brief . */
  BLAS_Cpp::Uplo	uplo() const {
    return uplo_;
  }
  /// Allow address to be taken of rvalue of this object
  DMatrixSliceTriEle* operator&() {
    return this;
  }
  /** \brief . */
  const DMatrixSliceTriEle* operator&() const {
    return this;
  }

private:	
  DMatrixSlice	gms_;
  BLAS_Cpp::Uplo	uplo_;
  // Not defined and not to be called
  DMatrixSliceTriEle();
  DMatrixSliceTriEle& operator=(const DMatrixSliceTriEle&);
};	// end class DMatrixSliceTriEle

inline
/// Return a triangular element-wise matrix
DMatrixSliceTriEle nonconst_tri_ele(DMatrixSlice gms, BLAS_Cpp::Uplo uplo)
{
  return DMatrixSliceTriEle(gms, uplo);
}

inline
/** \brief . */
const DMatrixSliceTriEle tri_ele(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
{
  return DMatrixSliceTriEle(gms, uplo);
}

// ////////////////////////////////////////////////////////////////////////////////
/** \brief . */
/* * Aggregate information for a triangular matrix (structure dependent) stored in a DMatrix.
  *
  * This is the type to be used as a rhs argument in linear algebra operations
  * that are structure specific like the BLAS operations.
  */
class DMatrixSliceTri {
public:
  /** \brief . */
  DMatrixSliceTri(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo, BLAS_Cpp::Diag diag)
    : gms_(const_cast<DMatrixSlice&>(gms)), uplo_(uplo), diag_(diag)
  {
    #ifdef LINALGPACK_CHECK_RHS_SIZES
      assert_gms_square(gms);
    #endif
  }
  /** \brief . */
  size_type rows() const {
    return gms_.rows();
  }
  /** \brief . */
  size_type cols() const {
    return gms_.cols();
  }
  /** \brief . */
  DMatrixSlice& gms() {
    return gms_;
  }
  /** \brief . */
  const DMatrixSlice& gms() const {
    return gms_;
  }
  /** \brief . */
  BLAS_Cpp::Uplo	uplo() const {
    return uplo_;
  }
  /** \brief . */
  BLAS_Cpp::Diag	diag() const {
    return diag_;
  }
  /// Allow address to be taken of rvalue of this object
  DMatrixSliceTri* operator&() {
    return this;
  }
  /** \brief . */
  const DMatrixSliceTri* operator&() const {
    return this;
  }

private:	
  DMatrixSlice	gms_;
  BLAS_Cpp::Uplo	uplo_;
  BLAS_Cpp::Diag	diag_;
  // not defined and not to be called
  DMatrixSliceTri();
  DMatrixSliceTri& operator=(const DMatrixSliceTri&);
};	// end class DMatrixSliceTri

inline
/// Return a triangular matrix
DMatrixSliceTri nonconst_tri(DMatrixSlice gms, BLAS_Cpp::Uplo uplo, BLAS_Cpp::Diag diag)
{
  return DMatrixSliceTri(gms, uplo, diag);
}

inline
/** \brief . */
const DMatrixSliceTri tri(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo, BLAS_Cpp::Diag diag)
{
  return DMatrixSliceTri(gms, uplo, diag);
}

// /////////////////////////////////////////////////////////////////////////////////////////
/** \brief . */
/* * Aggregate information for a symmetric matrix stored in a DMatrix.
  *
  * This is the type to be used as both lhs and rhs arguments in linear algebra operations.
  */
class DMatrixSliceSym {
public:
  /** \brief . */
  DMatrixSliceSym(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
    : gms_(const_cast<DMatrixSlice&>(gms)), uplo_(uplo)
  {
    #ifdef LINALGPACK_CHECK_RHS_SIZES
      assert_gms_square(gms);
    #endif
  }
  /** \brief . */
  size_type rows() const {
    return gms_.rows();
  }
  /** \brief . */
  size_type cols() const {
    return gms_.cols();
  }
  /** \brief . */
  DMatrixSlice& gms() {
    return gms_;
  }
  /** \brief . */
  const DMatrixSlice& gms() const {
    return gms_;
  }
  /** \brief . */
  BLAS_Cpp::Uplo	uplo() const {
    return uplo_;
  }
  /// Allow address to be taken of rvalue of this object
  DMatrixSliceSym* operator&() {
    return this;
  }
  const DMatrixSliceSym* operator&() const {
    return this;
  }

private:	
  DMatrixSlice	gms_;
  BLAS_Cpp::Uplo	uplo_;
  // not defined and not to be called
  DMatrixSliceSym();
  DMatrixSliceSym& operator=(const DMatrixSliceTri&);
};	// end class DMatrixSliceSym

inline
/// Return a symmetric matrix
DMatrixSliceSym nonconst_sym(DMatrixSlice gms, BLAS_Cpp::Uplo uplo)
{
  return DMatrixSliceSym(gms, uplo);
}

inline
/** \brief . */
const DMatrixSliceSym sym(const DMatrixSlice& gms, BLAS_Cpp::Uplo uplo)
{
  return DMatrixSliceSym(gms, uplo);
}

// @}

} // end namespace DenseLinAlgPack

#endif	// GEN_MATRIX_AS_TRI_SYM_H
