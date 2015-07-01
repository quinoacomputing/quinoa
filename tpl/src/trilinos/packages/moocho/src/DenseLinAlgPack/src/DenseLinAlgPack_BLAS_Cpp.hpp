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
// C++ overloads for BLAS kernals (element type removed from name and enum for operations)

#ifndef BLAS_CPP_OVERLOADS_DECLARATIONS_H
#define BLAS_CPP_OVERLOADS_DECLARATIONS_H

#include "Teuchos_F77_wrappers.h"
#include "BLAS_Cpp_Types.hpp"

// Overloaded BLAS wrappers.
// The naming convention is the Fortran BLAS name minus the type prefix.
namespace BLAS_Cpp {

typedef FortranTypes::f_int			f_int;
typedef FortranTypes::f_real		f_real;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;

/* * @name Option Arguments
 * These are emumerations that are used with the overloaded C++ BLAS declarations to replace the
 * error prone use of characters for specifying options
 * @memo enumerations (enum)
 */
// @{

/** \brief . */
const char SideChar[]	= {'L'	, 'R'			};
/** \brief . */
const char TransChar[]	= {'N'	, 'T'	, 'C'	};
/** \brief . */
const char UploChar[]	= {'U'	, 'L'			};
/** \brief . */
const char DiagChar[]	= {'U'	, 'N'			};

// @}

/* * @name C++ BLAS Function Declarations
 * These are overloaded C++ functions that have removed the element type from the name
 * of the BLAS functions and use enumerations for the options arguments.
 */
// @{	

// ///////////////////////////////////////////////////////////////////////////////////////////
/* * @name Level 1 BLAS (vector-vector operations) */
// @{	

/* * @name Generate plane rotation */
// @{

/** \brief . */
void rotg( f_dbl_prec* a, f_dbl_prec* b, f_dbl_prec* c, f_dbl_prec* s );

// @}	 	
 
/* * @name Apply plane rotation */
// @{

/** \brief . */
void rot(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY
     , const f_dbl_prec& C, const f_dbl_prec& S);
// @}

/* * @name  Interchange vectors */
// @{

/** \brief . */
void swap(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY);

// @}

/* * @name  DVector scaling */
// @{

/// 
void scal(const f_int& N, const f_dbl_prec& ALPHA, f_dbl_prec* X, const f_int& INCX);

// @}

/* * @name DVector copy */
// @{

/// 
void copy(const f_int& N, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY);

// @}

/* * @name  y = a*x + y */
// @{

/** \brief . */
void axpy(const f_int& N, const f_dbl_prec& A, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y
  , const f_int& INCY);
  
// @}

/* * @name  Dot product */
// @{

/** \brief . */
f_dbl_prec dot(const f_int& N, const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec* Y, const f_int& INCY);

// @}

/* * @name  2-Norm */
// @{

/** \brief . */
f_dbl_prec nrm2(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// @}

/* * @name  1-Norm */
// @{

/** \brief . */
f_dbl_prec asum(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// @}

/* * @name  Inifinity-Norm */
// @{

/** \brief . */
f_dbl_prec iamax(const f_int& N, const f_dbl_prec* X, const f_int& INCX);
// @}

//		end Level-1 BLAS
// @}	

// /////////////////////////////////////////////////
/* * @name Level-2 BLAS (matrix-vector operations) */
// @{	

/* * @name General rectangular matrix-vector products */
// @{

/** \brief . */
void gemv(Transp transa, f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);
       
// @}	

/* * @name General band matrix-vector products */
// @{

/** \brief . */
void gbmv(Transp transa, f_int m, f_int n, f_int kl, f_int ku, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);
         
// @}

/* * @name Hermitian matrix-vector products */
// @{

         
// @}

/* * @name Hermitian band matrix-vector products */
// @{

// @}

/* * @name Hermitian packed matrix-vector products */
// @{


// @}

/* * @name Symmetric matrix-vector products */
// @{

/** \brief . */
void symv(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);
      
// @}

/* * @name Symmetric band matrix-vector products */
// @{

/** \brief . */
void sbmv(Uplo uplo, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);

// @}

/* * @name Symmetric packed matrix-vector products */
// @{

/** \brief . */
void spmv(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* pap
  , const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);

// @}

/* * @name Triangular matrix-vector products */
// @{

/** \brief . */
void trmv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular band matrix-vector products */
// @{

/** \brief . */
void tbmv(Uplo uplo, Transp trans, Diag diag, f_int n, f_int k, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular packed matrix-vector products */
// @{

/** \brief . */
void tpmv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pap
  , f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular equation solve */
// @{

/** \brief . */
void trsv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular band equation solve */
// @{

/** \brief . */
void tbsv(Uplo uplo, Transp trans, Diag diag, f_int n, f_int k, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular packed equation solve */
// @{

/** \brief . */
void tpsv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pap
  , f_dbl_prec* px, f_int incx);

// @}

/* * @name General rank-1 update */
// @{

/** \brief . */
void ger(f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pa, f_int lda);

// @}

/* * @name Hermitian rank-1 update */
// @{

// @}

/* * @name Hermitian packed rank-1 update */
// @{

// @}

/* * @name Hermitian rank-2 update */
// @{

// @}

/* * @name Hermitian packed rank-2 update */
// @{

// @}

/* * @name Symmetric rank-1 update */
// @{

/** \brief . */
void syr(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, f_dbl_prec* pa, f_int lda);

// @}

/* * @name Symmetric packed rank-1 update */
// @{

/** \brief . */
void spr(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, f_dbl_prec* pap);

// @}

/* * @name Symmetric rank-2 update */
// @{

/** \brief . */
void syr2(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pa, f_int lda);

// @}

/* * @name Symmetric packed rank-2 update */
// @{

/** \brief . */
void spr2(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pap);

// @}

//		end Level 2 BLAS
// @}
  
// /////////////////////////////////////////
/* * @name Level 3 BLAS (matrix-matrix operations) */
// @{	

/* * @name General rectangular matrix-matrix product */
// @{

/** \brief . */
void gemm(Transp transa, Transp transb, f_int m, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc);

// @}

/* * @name Symmetric matrix-matrix product */
// @{

/** \brief . */
void symm(Side side, Uplo uplo, f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc);

// @}

/* * @name Hermitian matrix-matrix product */
// @{

// @}

/* * @name Symmetric rank-k update */
// @{

/** \brief . */
void syrk(Uplo uplo, Transp trans, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc);

// @}

/* * @name Hermitian rank-k update */
// @{

// @}

/* * @name Symmetric rank-2k update */
// @{

/** \brief . */
void syr2k(Uplo uplo, Transp trans, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc);

// @}

/* * @name Hermitian rank-2k update */
// @{

// @}

/* * @name Triangular matrix-matrix product */
// @{

/** \brief . */
void trmm(Side side, Uplo uplo, Transp transa, Diag diag, f_int m, f_int n, f_dbl_prec alpha
  , const f_dbl_prec* pa, f_int lda, f_dbl_prec* pb, f_int ldb);

// @}

/* * @name Solution of triangular system */
// @{

/** \brief . */
void trsm(Side side, Uplo uplo, Transp transa, Diag diag, f_int m, f_int n, f_dbl_prec alpha
  , const f_dbl_prec* pa, f_int lda, f_dbl_prec* pb, f_int ldb);
          
// @}

//		end Level 3 BLAS
// @}	

//		end overloaded functions
// @}

}	// end namespace BLAS_Cpp

#endif // BLAS_CPP_OVERLOADS_DECLARATIONS_H
