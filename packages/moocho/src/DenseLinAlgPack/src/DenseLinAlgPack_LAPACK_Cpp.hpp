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

#ifndef LAPACK_CPP_H
#define LAPACK_CPP_H

#include "DenseLinAlgPack_LAPACK_C_Decl.hpp"
#include "DenseLinAlgPack_BLAS_Cpp.hpp"

// Cpp Declarations for calling LAPACK functions that
// use function overloading to remove the floating point
// typeds from the names and use enumerations to replace
// the char arguments.

namespace LAPACK_Cpp {

typedef FortranTypes::f_int			f_int;
typedef FortranTypes::f_real		f_real;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;
typedef FortranTypes::f_logical		f_logical;

// xPOTRF

inline
/** \brief . */
void potrf(  BLAS_Cpp::Uplo uplo
  , const f_int& n, f_dbl_prec* A, const f_int& lda
  , f_int* info )
{
  LAPACK_C_Decl::dpotrf( BLAS_Cpp::UploChar[uplo]
    ,n ,A ,lda ,info );
}

// xGEQRF

inline
/** \brief . */
void geqrf( const f_int& m
  , const f_int& n, f_dbl_prec* A, const f_int& lda
  , f_dbl_prec* tau, f_dbl_prec* work
  , const f_int& lwork, f_int* info  )
{
  LAPACK_C_Decl::dgeqrf(m,n,A,lda,tau,work,lwork,info);
}

// xORMQR

inline
/** \brief . */
void ormqr( BLAS_Cpp::Side side, BLAS_Cpp::Transp trans
  , const f_int& m, const f_int& n
  , const f_int& k, const f_dbl_prec* A, const f_int& lda
  , const f_dbl_prec* tau, f_dbl_prec* C, const f_int& ldc
  , f_dbl_prec* work, const f_int& lwork, f_int* info )
{
  LAPACK_C_Decl::dormqr( BLAS_Cpp::SideChar[side]
    , BLAS_Cpp::TransChar[trans], m, n, k, A, lda
    , tau, C, ldc, work, lwork, info );
  
}

// xSYTRF

inline
//
void sytrf( BLAS_Cpp::Uplo uplo, const f_int& n, f_dbl_prec A[]
  , const f_int& lda, f_int ipiv[], f_dbl_prec work[], const f_int& lwork
  , f_int* info )
{
  LAPACK_C_Decl::dsytrf( BLAS_Cpp::UploChar[uplo]
    , n, A, lda, ipiv, work, lwork, info );
}

// xSYTRS

inline
/** \brief . */
void sytrs( BLAS_Cpp::Uplo uplo
  , const f_int& n, const f_int& nrhs, const f_dbl_prec A[]
  , const f_int& lda, const f_int ipiv[], f_dbl_prec B[]
  , const f_int& ldb, f_int* info )
{
  LAPACK_C_Decl::dsytrs( BLAS_Cpp::UploChar[uplo]
    , n, nrhs, A, lda, ipiv, B, ldb, info );
}

// xGETRF

inline
//
void getrf(
  const f_int& m, const f_int& n, f_dbl_prec A[]
  ,const f_int& lda, f_int ipiv[], f_int* info )
{
  LAPACK_C_Decl::dgetrf( m, n, A, lda, ipiv, info );
}

// xGETRS

inline
/** \brief . */
void getrs(
  BLAS_Cpp::Transp trans
  ,const f_int& n, const f_int& nrhs, const f_dbl_prec A[]
  , const f_int& lda, const f_int ipiv[], f_dbl_prec B[]
  , const f_int& ldb, f_int* info )
{
  LAPACK_C_Decl::dgetrs(
    BLAS_Cpp::TransChar[trans], n, nrhs, A, lda, ipiv, B, ldb, info
    );
}

}	// end namespace LAPACK_Cpp

#endif	// LAPACK_CPP_H
