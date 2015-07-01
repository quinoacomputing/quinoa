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

#include "DenseLinAlgPack_BLAS_Cpp.hpp"

// /////////////////////////////////////
// Fortran function declarations.

namespace BLAS_C_Decl {

typedef FortranTypes::f_int  		f_int;
typedef FortranTypes::f_real		f_real;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;

// function declarations 
// (don't use these directly, use the later overloaded wrappers in namespace BLAS)
extern "C" {

// ////////////////////////////////////////
// Level 1 BLAS (vector-vector operations)

// Generate plane rotation
FORTRAN_FUNC_DECL_UL(void,DROTG,drotg)(f_dbl_prec* A, f_dbl_prec* B, f_dbl_prec* C, f_dbl_prec* S);
 
// Apply plane rotation
FORTRAN_FUNC_DECL_UL(void,DROT,drot)(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY,
      const f_dbl_prec& C, const f_dbl_prec& S);

// Interchange vectors
FORTRAN_FUNC_DECL_UL(void,DSWAP,dswap)(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY);

// DVector scaling
FORTRAN_FUNC_DECL_UL(void,DSCAL,dscal)(const f_int& N, const f_dbl_prec& ALPHA, f_dbl_prec* X, const f_int& INCX);

// DVector copy 
FORTRAN_FUNC_DECL_UL(void,DCOPY,dcopy)(const f_int& N, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY);

// y = a*x + y
FORTRAN_FUNC_DECL_UL(void,DAXPY,daxpy)(const f_int& N, const f_dbl_prec& A, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y,
       const f_int& INCY);

// Dot product
FORTRAN_FUNC_DECL_UL(f_dbl_prec,DDOT,ddot)(const f_int& N, const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec* Y, const f_int& INCY);
FORTRAN_FUNC_DECL_UL(f_dbl_prec,DSDOT,dsdot)(const f_int& N, const f_real* X, const f_int& INCX, const f_real* Y, const f_int& INCY);

// 2-Norm
FORTRAN_FUNC_DECL_UL(f_dbl_prec,DNRM2,dnrm2)(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// Sum of magnitudes
FORTRAN_FUNC_DECL_UL(f_dbl_prec,DASUM,dasum)(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// Largest component of vector
FORTRAN_FUNC_DECL_UL(f_dbl_prec,IDAMAX,idamax)(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// ////////////////////////////////////////
// Level 2 BLAS (matrix-vector operations)

// General rectangular matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DGEMV,dgemv)(FORTRAN_CONST_CHAR_1_ARG(TRANSA)
  , const f_int& M, const f_int& N, const f_dbl_prec& ALPHA
  , const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* X, const f_int& INCX
  , const f_dbl_prec& BETA, f_dbl_prec* Y, const f_int& INCY
  );

// General band matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DGBMV,dgbmv)(FORTRAN_CONST_CHAR_1_ARG(TRANSA)
  , const f_int& M, const f_int& N, const f_int& KL, const f_int& KU
  , const f_dbl_prec& ALPHA,	const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* X
  , const f_int& INCX, const f_dbl_prec& BETA, f_dbl_prec* Y, const f_int& INCY
  );

// Hermitian matrix-vector products

// Hermitian band matrix-vector products

// Hermitian packed matrix-vector products

// Symmetric matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DSYMV,dsymv)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , const f_int& N
  , const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA
  , const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec& BETA
  , f_dbl_prec* Y, const f_int& INCY
  );

// Symmetric band matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DSBMV,dsbmv)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , const f_int& N, const f_int& K
  , const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA
  , const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec& BETA
  , f_dbl_prec* Y, const f_int& INCY
  );

// Symmetric packed matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DSPMV,dspmv)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , const f_int& N
  , const f_dbl_prec& ALPHA, const f_dbl_prec* AP
  , const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec& BETA
  , f_dbl_prec* Y, const f_int& INCY
  );

// Triangular matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DTRMV,dtrmv)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , FORTRAN_CONST_CHAR_1_ARG(TRANS), FORTRAN_CONST_CHAR_1_ARG(DIAG), const f_int& N
  , const f_dbl_prec* A, const f_int& LDA, f_dbl_prec* X, const f_int& INCX);

// Triangular band matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DTBMV,dtbmv)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , FORTRAN_CONST_CHAR_1_ARG(TRANS), FORTRAN_CONST_CHAR_1_ARG(DIAG)
  , const f_int& N, const f_int& K
  , const f_dbl_prec* A, const f_int& LDA, f_dbl_prec* X, const f_int& INCX);

// Triangular packed matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DTPMV,dtpmv)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , FORTRAN_CONST_CHAR_1_ARG(TRANS), FORTRAN_CONST_CHAR_1_ARG(DIAG), const f_int& N
  , const f_dbl_prec* AP, f_dbl_prec* X, const f_int& INCX);

// Triangular equation solve
FORTRAN_FUNC_DECL_UL(void,DTRSV,dtrsv)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , FORTRAN_CONST_CHAR_1_ARG(TRANS), FORTRAN_CONST_CHAR_1_ARG(DIAG), const f_int& N
  , const f_dbl_prec* A, const f_int& LDA, f_dbl_prec* X, const f_int& INCX);

// Triangular band equation solve
FORTRAN_FUNC_DECL_UL(void,DTBSV,dtbsv)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , FORTRAN_CONST_CHAR_1_ARG(TRANS), FORTRAN_CONST_CHAR_1_ARG(DIAG), const f_int& N
  , const f_int& K
  , const f_dbl_prec* A, const f_int& LDA, f_dbl_prec* X, const f_int& INCX);

// Triangular packed equation solve
FORTRAN_FUNC_DECL_UL(void,DTPSV,dtpsv)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , FORTRAN_CONST_CHAR_1_ARG(TRANS), FORTRAN_CONST_CHAR_1_ARG(DIAG), const f_int& N
  , const f_dbl_prec* AP, f_dbl_prec* X, const f_int& INCX);

// General rank-1 update
FORTRAN_FUNC_DECL_UL(void,DGER,dger)(const f_int& M, const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX,
      const f_dbl_prec* Y, const f_int& INCY, f_dbl_prec* A, const f_int& LDA);

// Hermitian rank-1 update

// Hermitian packed rank-1 update

// Hermitian rank-2 update

// Hermitian packed rank-2 update

// Symmetric rank-1 update
FORTRAN_FUNC_DECL_UL(void,DSYR,dsyr)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX
  , f_dbl_prec* A, const f_int& LDA);

// Symmetric packed rank-1 update
FORTRAN_FUNC_DECL_UL(void,DSPR,dspr)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX
  , f_dbl_prec* AP);

// Symmetric rank-2 update
FORTRAN_FUNC_DECL_UL(void,DSYR2,dsyr2)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX
  , const f_dbl_prec* Y, const f_int& INCY, f_dbl_prec* A, const f_int& LDA);

// Symmetric packed rank-2 update
FORTRAN_FUNC_DECL_UL(void,DSPR2,dspr2)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX
  , const f_dbl_prec* Y, const f_int& INCY, f_dbl_prec* AP);

// ////////////////////////////////////////
// Level 3 BLAS (matrix-matrix operations)

// General rectangular matrix-matrix product
FORTRAN_FUNC_DECL_UL(void,DGEMM,dgemm)(FORTRAN_CONST_CHAR_1_ARG(TRANSA)
  , FORTRAN_CONST_CHAR_1_ARG(TRANSB), const f_int& M, const f_int& N, const f_int& K
  , const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* B
  , const f_int& LDB, const f_dbl_prec& BETA, f_dbl_prec* C, const f_int& LDC);

// Symmetric matrix-matrix product
FORTRAN_FUNC_DECL_UL(void,DSYMM,dsymm)(FORTRAN_CONST_CHAR_1_ARG(SIDE)
  , FORTRAN_CONST_CHAR_1_ARG(UPLO), const f_int& M, const f_int& N
  , const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* B
  , const f_int& LDB, const f_dbl_prec& BETA, f_dbl_prec* C, const f_int& LDC);

// Hermitian matrix-matrix product

// Symmetric rank-k update
FORTRAN_FUNC_DECL_UL(void,DSYRK,dsyrk)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , FORTRAN_CONST_CHAR_1_ARG(TRANS), const f_int& N, const f_int& K
  , const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA
  , const f_dbl_prec& BETA, f_dbl_prec* C, const f_int& LDC);

// Hermitian rank-k update

// Symmetric rank-2k update
FORTRAN_FUNC_DECL_UL(void,DSYR2K,dsyr2k)(FORTRAN_CONST_CHAR_1_ARG(UPLO)
  , FORTRAN_CONST_CHAR_1_ARG(TRANS), const f_int& N, const f_int& K
  , const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* B
  , const f_int& LDB, const f_dbl_prec& BETA, f_dbl_prec* C, const f_int& LDC);

// Hermitian rank-2k update

// Triangular matrix-matrix product
FORTRAN_FUNC_DECL_UL(void,DTRMM,dtrmm)(FORTRAN_CONST_CHAR_1_ARG(SIDE)
  , FORTRAN_CONST_CHAR_1_ARG(UPLO), FORTRAN_CONST_CHAR_1_ARG(TRANSA)
  , FORTRAN_CONST_CHAR_1_ARG(DIAG)
  , const f_int& M, const f_int& N,	const f_dbl_prec& ALPHA, const f_dbl_prec* A
  , const f_int& LDA, f_dbl_prec* B, const f_int& LDB);

// Solution of triangular system 
FORTRAN_FUNC_DECL_UL(void,DTRSM,dtrsm)(FORTRAN_CONST_CHAR_1_ARG(SIDE)
  , FORTRAN_CONST_CHAR_1_ARG(UPLO), FORTRAN_CONST_CHAR_1_ARG(TRANSA)
  , FORTRAN_CONST_CHAR_1_ARG(DIAG)
  , const f_int& M, const f_int& N,	const f_dbl_prec& ALPHA
  , const f_dbl_prec* A, const f_int& LDA, f_dbl_prec* B, const f_int& LDB);

} // end extern "C"

} // end namespace BLAS_C_Decl

// ////////////////////////////////////////////////////////
// C++ BLAS Function Declarations

// Level 1 BLAS (vector-vector operations)

void BLAS_Cpp::rotg(f_dbl_prec* a, f_dbl_prec* b, f_dbl_prec* c, f_dbl_prec* s) {
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DROTG,drotg)(a,b,c,s);
}
 
// Apply plane rotation

void BLAS_Cpp::rot(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY
  , const f_dbl_prec& C, const f_dbl_prec& S)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DROT,drot)(N, X, INCX, Y, INCY, C, S);
}

// Interchange vectors

void BLAS_Cpp::swap(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSWAP,dswap)(N, X, INCX, Y, INCY);
}		 	

// DVector scaling

void BLAS_Cpp::scal(const f_int& N, const f_dbl_prec& ALPHA, f_dbl_prec* X, const f_int& INCX)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSCAL,dscal)(N, ALPHA, X, INCX);
}

// DVector copy

void BLAS_Cpp::copy(const f_int& N, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DCOPY,dcopy)(N, X, INCX, Y, INCY);
}

// y = a*x + y

void BLAS_Cpp::axpy(const f_int& N, const f_dbl_prec& A, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y
  , const f_int& INCY)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DAXPY,daxpy)(N, A, X, INCX, Y, INCY);
}

// Dot product

BLAS_Cpp::f_dbl_prec BLAS_Cpp::dot(const f_int& N, const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec* Y, const f_int& INCY)
{
  return BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DDOT,ddot)(N, X, INCX, Y, INCY);
}

// 2-Norm

BLAS_Cpp::f_dbl_prec BLAS_Cpp::nrm2(const f_int& N, const f_dbl_prec* X, const f_int& INCX)
{
  return BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DNRM2,dnrm2)(N, X, INCX);
}

// 1-Norm

BLAS_Cpp::f_dbl_prec BLAS_Cpp::asum(const f_int& N, const f_dbl_prec* X, const f_int& INCX)
{
  return BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DASUM,dasum)(N, X, INCX);
}

// Inifinity-Norm

BLAS_Cpp::f_dbl_prec BLAS_Cpp::iamax(const f_int& N, const f_dbl_prec* X, const f_int& INCX)
{
  return BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(IDAMAX,idamax)(N, X, INCX);
}

// Level-2 BLAS (matrix-vector operations)

// General rectangular matrix-vector products

void BLAS_Cpp::gemv(Transp transa, f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DGEMV,dgemv)(FORTRAN_CHAR_1_ARG_CALL(TransChar[transa]), m, n, alpha, pa, lda, x, incx, beta, py, incy);
}

// General band matrix-vector products

void BLAS_Cpp::gbmv(Transp transa, f_int m, f_int n, f_int kl, f_int ku, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DGBMV,dgbmv)(FORTRAN_CHAR_1_ARG_CALL(TransChar[transa]), m, n, kl, ku, alpha, pa, lda, x, incx, beta, py, incy);
}
         
// Hermitian matrix-vector products

// Hermitian band matrix-vector products

// Hermitian packed matrix-vector products

// Symmetric matrix-vector products

void BLAS_Cpp::symv(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYMV,dsymv)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), n, alpha, pa, lda, x, incx, beta, py, incy);
}
      
// Symmetric band matrix-vector products

void BLAS_Cpp::sbmv(Uplo uplo, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSBMV,dsbmv)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), n, k, alpha, pa, lda, x, incx, beta, py, incy);
}

// Symmetric packed matrix-vector products

void BLAS_Cpp::spmv(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* pap
  , const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSPMV,dspmv)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), n, alpha, pap, x, incx, beta, py, incy);
}

// Triangular matrix-vector products

void BLAS_Cpp::trmv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec* px, f_int incx)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTRMV,dtrmv)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo])
    ,FORTRAN_CHAR_1_ARG_CALL(TransChar[trans]),FORTRAN_CHAR_1_ARG_CALL(DiagChar[diag])
    ,n, pa, lda, px, incx);
}

// Triangular band matrix-vector products

void BLAS_Cpp::tbmv(Uplo uplo, Transp trans, Diag diag, f_int n, f_int k, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec* px, f_int incx)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTBMV,dtbmv)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo])
    ,FORTRAN_CHAR_1_ARG_CALL(TransChar[trans]), FORTRAN_CHAR_1_ARG_CALL(DiagChar[diag])
    ,n, k, pa, lda, px, incx);
}

// Triangular packed matrix-vector products

void BLAS_Cpp::tpmv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pap
  , f_dbl_prec* px, f_int incx)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTPMV,dtpmv)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo])
    ,FORTRAN_CHAR_1_ARG_CALL(TransChar[trans]), FORTRAN_CHAR_1_ARG_CALL(DiagChar[diag])
    ,n, pap, px, incx);
}

// Triangular equation solve

void BLAS_Cpp::trsv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec* px, f_int incx)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTRSV,dtrsv)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo])
    ,FORTRAN_CHAR_1_ARG_CALL(TransChar[trans]), FORTRAN_CHAR_1_ARG_CALL(DiagChar[diag])
    ,n, pa, lda, px, incx);
}

// Triangular band equation solve

void BLAS_Cpp::tbsv(Uplo uplo, Transp trans, Diag diag, f_int n, f_int k, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec* px, f_int incx)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTBSV,dtbsv)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo])
    ,FORTRAN_CHAR_1_ARG_CALL(TransChar[trans]), FORTRAN_CHAR_1_ARG_CALL(DiagChar[diag])
    ,n, k, pa, lda, px, incx);
}

// Triangular packed equation solve

void BLAS_Cpp::tpsv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pap
  , f_dbl_prec* px, f_int incx)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTPSV,dtpsv)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo])
    ,FORTRAN_CHAR_1_ARG_CALL(TransChar[trans]), FORTRAN_CHAR_1_ARG_CALL(DiagChar[diag])
    ,n, pap, px, incx);
}

// General rank-1 update

void BLAS_Cpp::ger(f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pa, f_int lda)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DGER,dger)(m, n, alpha, px, incx, py, incy, pa, lda);
}

// Hermitian rank-1 update

// Hermitian packed rank-1 update

// Hermitian rank-2 update

// Hermitian packed rank-2 update

// Symmetric rank-1 update

void BLAS_Cpp::syr(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, f_dbl_prec* pa, f_int lda)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYR,dsyr)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), n, alpha, px, incx, pa, lda);
}

// Symmetric packed rank-1 update

void BLAS_Cpp::spr(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, f_dbl_prec* pap)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSPR,dspr)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), n, alpha, px, incx, pap);
}

// Symmetric rank-2 update

void BLAS_Cpp::syr2(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pa, f_int lda)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYR2,dsyr2)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), n, alpha, px, incx, py, incy, pa, lda);
}

// Symmetric packed rank-2 update

void BLAS_Cpp::spr2(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
  , f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pap)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSPR2,dspr2)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), n, alpha, px, incx, py, incy, pap);
}

// Level 3 BLAS (matrix-matrix operations)	

// General rectangular matrix-matrix product

void BLAS_Cpp::gemm(Transp transa, Transp transb, f_int m, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DGEMM,dgemm)(FORTRAN_CHAR_1_ARG_CALL(TransChar[transa])
    ,FORTRAN_CHAR_1_ARG_CALL(TransChar[transb]), m, n, k, alpha, pa, lda, pb, ldb
    , beta, pc, ldc);
}

// Symmetric matrix-matrix product

void BLAS_Cpp::symm(Side side, Uplo uplo, f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYMM,dsymm)(FORTRAN_CHAR_1_ARG_CALL(SideChar[side]), FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}

// Hermitian matrix-matrix product

// Symmetric rank-k update

void BLAS_Cpp::syrk(Uplo uplo, Transp trans, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYRK,dsyrk)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo])
    ,FORTRAN_CHAR_1_ARG_CALL(TransChar[trans]), n, k, alpha, pa, lda, beta, pc, ldc);
}

// Hermitian rank-k update

// Symmetric rank-2k update

void BLAS_Cpp::syr2k(Uplo uplo, Transp trans, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
  , f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYR2K,dsyr2k)(FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo])
    ,FORTRAN_CHAR_1_ARG_CALL(TransChar[trans]), n, k, alpha, pa, lda, pb, ldb
    ,beta, pc, ldc);
}

// Hermitian rank-2k update

// Triangular matrix-matrix product

void BLAS_Cpp::trmm(Side side, Uplo uplo, Transp transa, Diag diag, f_int m, f_int n, f_dbl_prec alpha
  , const f_dbl_prec* pa, f_int lda, f_dbl_prec* pb, f_int ldb)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTRMM,dtrmm)(FORTRAN_CHAR_1_ARG_CALL(SideChar[side])
    ,FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), FORTRAN_CHAR_1_ARG_CALL(TransChar[transa])
    ,FORTRAN_CHAR_1_ARG_CALL(DiagChar[diag]), m, n, alpha, pa, lda, pb, ldb);
}

// Solution of triangular system

void BLAS_Cpp::trsm(Side side, Uplo uplo, Transp transa, Diag diag, f_int m, f_int n, f_dbl_prec alpha
  , const f_dbl_prec* pa, f_int lda, f_dbl_prec* pb, f_int ldb)
{
  BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTRSM,dtrsm)(FORTRAN_CHAR_1_ARG_CALL(SideChar[side])
    ,FORTRAN_CHAR_1_ARG_CALL(UploChar[uplo]), FORTRAN_CHAR_1_ARG_CALL(TransChar[transa])
    ,FORTRAN_CHAR_1_ARG_CALL(DiagChar[diag]), m, n, alpha, pa, lda, pb, ldb);
}
