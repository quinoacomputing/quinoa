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

#include "DenseLinAlgPack_LAPACK_C_Decl.hpp"

// DPOTRF

extern "C" {
  FORTRAN_FUNC_DECL_UL( void, DPOTRF, dpotrf ) ( FORTRAN_CONST_CHAR_1_ARG(UPLO)
    , const LAPACK_C_Decl::f_int& N, LAPACK_C_Decl::f_dbl_prec* A
    , const LAPACK_C_Decl::f_int& LDA, LAPACK_C_Decl::f_int* INFO  );
}

void LAPACK_C_Decl::dpotrf(  const f_char& UPLO
  , const f_int& N, f_dbl_prec* A, const f_int& LDA
  , f_int* INFO  )
{
  ::FORTRAN_FUNC_CALL_UL(DPOTRF,dpotrf) (FORTRAN_CHAR_1_ARG_CALL(UPLO),N,A,LDA,INFO);
}

// DGEQRF

extern "C" {
  FORTRAN_FUNC_DECL_UL( void, DGEQRF, dgeqrf ) (
    const LAPACK_C_Decl::f_int& M, const LAPACK_C_Decl::f_int& N
    , LAPACK_C_Decl::f_dbl_prec* A, const LAPACK_C_Decl::f_int& LDA
    , LAPACK_C_Decl::f_dbl_prec* TAU, LAPACK_C_Decl::f_dbl_prec* WORK
    , const LAPACK_C_Decl::f_int& LWORK, LAPACK_C_Decl::f_int* INFO  );
}

void LAPACK_C_Decl::dgeqrf( const f_int& M
  , const f_int& N, f_dbl_prec* A, const f_int& LDA
  , f_dbl_prec* TAU, f_dbl_prec* WORK
  , const f_int& LWORK, f_int* INFO  )
{
  ::FORTRAN_FUNC_CALL_UL(DGEQRF,dgeqrf) (M,N,A,LDA,TAU,WORK,LWORK,INFO);
}

// DORMRQ

extern "C" {
  FORTRAN_FUNC_DECL_UL( void, DORMQR, dormqr ) ( FORTRAN_CONST_CHAR_1_ARG(SIDE)
    , FORTRAN_CONST_CHAR_1_ARG(TRANS)
    , const LAPACK_C_Decl::f_int& M, const LAPACK_C_Decl::f_int& N
    , const LAPACK_C_Decl::f_int& K, const LAPACK_C_Decl::f_dbl_prec* A
    , const LAPACK_C_Decl::f_int& LDA, const LAPACK_C_Decl::f_dbl_prec* TAU
    , LAPACK_C_Decl::f_dbl_prec* C, const LAPACK_C_Decl::f_int& LDC
    , LAPACK_C_Decl::f_dbl_prec* WORK, const LAPACK_C_Decl::f_int& LWORK
    , LAPACK_C_Decl::f_int* INFO );
}

void LAPACK_C_Decl::dormqr( const f_char& SIDE
  , const f_char& TRANS, const f_int& M, const f_int& N
  , const f_int& K, const f_dbl_prec* A, const f_int& LDA
  , const f_dbl_prec* TAU, f_dbl_prec* C, const f_int& LDC
  , f_dbl_prec* WORK, const f_int& LWORK, f_int* INFO )
{
  ::FORTRAN_FUNC_CALL_UL(DORMQR,dormqr)(FORTRAN_CHAR_1_ARG_CALL(SIDE)
    ,FORTRAN_CHAR_1_ARG_CALL(TRANS),M,N,K,A,LDA,TAU,C,LDC,WORK,LWORK,INFO);
}

// DSYTRF

extern "C" {
  FORTRAN_FUNC_DECL_UL( void, DSYTRF, dsytrf ) ( FORTRAN_CONST_CHAR_1_ARG(UPLO)
    , const LAPACK_C_Decl::f_int& N, LAPACK_C_Decl::f_dbl_prec A[]
    , const LAPACK_C_Decl::f_int& LDA
    , LAPACK_C_Decl::f_int IPIV[], LAPACK_C_Decl::f_dbl_prec WORK[]
    , const LAPACK_C_Decl::f_int& LWORK
    , LAPACK_C_Decl::f_int* INFO );
}

void LAPACK_C_Decl::dsytrf( const f_char& UPLO
  , const f_int& N, f_dbl_prec A[], const f_int& LDA
  , f_int IPIV[], f_dbl_prec WORK[], const f_int& LWORK
  , f_int* INFO )
{
  ::FORTRAN_FUNC_CALL_UL(DSYTRF,dsytrf)(FORTRAN_CHAR_1_ARG_CALL(UPLO)
    ,N,A,LDA,IPIV,WORK,LWORK,INFO);
}

// DSYTRS

extern "C" {
  FORTRAN_FUNC_DECL_UL( void, DSYTRS, dsytrs ) ( FORTRAN_CONST_CHAR_1_ARG(UPLO)
    , const LAPACK_C_Decl::f_int& N, const LAPACK_C_Decl::f_int& NRHS
    , const LAPACK_C_Decl::f_dbl_prec A[]
    , const LAPACK_C_Decl::f_int& LDA, const LAPACK_C_Decl::f_int IPIV[]
    , LAPACK_C_Decl::f_dbl_prec B[]
    , const LAPACK_C_Decl::f_int& LDB, LAPACK_C_Decl::f_int* INFO );
}

void LAPACK_C_Decl::dsytrs( const f_char& UPLO
  , const f_int& N, const f_int& NRHS, const f_dbl_prec A[]
  , const f_int& LDA, const f_int IPIV[], f_dbl_prec B[]
  , const f_int& LDB, f_int* INFO )
{
  ::FORTRAN_FUNC_CALL_UL(DSYTRS,dsytrs)(FORTRAN_CHAR_1_ARG_CALL(UPLO)
    ,N,NRHS,A,LDA,IPIV,B,LDB,INFO);
}

// DGETRF

extern "C" {
  FORTRAN_FUNC_DECL_UL( void, DGETRF, dgetrf ) (
    const LAPACK_C_Decl::f_int& M ,const LAPACK_C_Decl::f_int& N
    ,LAPACK_C_Decl::f_dbl_prec A[], const LAPACK_C_Decl::f_int& LDA
    ,LAPACK_C_Decl::f_int IPIV[], LAPACK_C_Decl::f_int* INFO
    );
}

void LAPACK_C_Decl::dgetrf(
  const f_int& M, const f_int& N, f_dbl_prec A[], const f_int& LDA
  ,f_int IPIV[], f_int* INFO
  )
{
  ::FORTRAN_FUNC_CALL_UL(DGETRF,dgetrf)(M,N,A,LDA,IPIV,INFO);
}

// DGETRS

extern "C" {
  FORTRAN_FUNC_DECL_UL( void, DGETRS, dgetrs ) (
    FORTRAN_CONST_CHAR_1_ARG(TRANS)
    ,const LAPACK_C_Decl::f_int& N, const LAPACK_C_Decl::f_int& NRHS
    ,const LAPACK_C_Decl::f_dbl_prec* A, const LAPACK_C_Decl::f_int& LDA
    ,const LAPACK_C_Decl::f_int IPIV[]
    ,LAPACK_C_Decl::f_dbl_prec* B, const LAPACK_C_Decl::f_int& LDB
    ,LAPACK_C_Decl::f_int* INFO
    );
}

void LAPACK_C_Decl::dgetrs(
  const f_char& TRANS
  ,const f_int& N, const f_int& NRHS, const f_dbl_prec A[]
  ,const f_int& LDA, const f_int IPIV[], f_dbl_prec B[]
  ,const f_int& LDB, f_int* INFO
  )
{
  ::FORTRAN_FUNC_CALL_UL(DGETRS,dgetrs)(
    FORTRAN_CHAR_1_ARG_CALL(TRANS),N,NRHS,A,LDA,IPIV,B,LDB,INFO
    );
}
