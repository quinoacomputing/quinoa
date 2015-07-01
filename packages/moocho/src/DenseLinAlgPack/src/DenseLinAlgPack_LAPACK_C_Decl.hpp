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

#ifndef LAPACK_C_DECL_H
#define LAPACK_C_DECL_H

#include "Teuchos_F77_wrappers.h"

// C Declarations for calling LAPACK functions.

namespace LAPACK_C_Decl {

typedef FortranTypes::f_int			f_int;
typedef FortranTypes::f_real		f_real;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;
typedef FortranTypes::f_logical		f_logical;
typedef FortranTypes::f_char		f_char;

// DPOTRF

void dpotrf(  const f_char& UPLO
  , const f_int& N, f_dbl_prec* A, const f_int& LDA
  , f_int* INFO  );

// DGEQRF

void dgeqrf( const f_int& M
  , const f_int& N, f_dbl_prec* A, const f_int& LDA
  , f_dbl_prec* TAU, f_dbl_prec* WORK
  , const f_int& LWORK, f_int* INFO  );

// DORMRQ

void dormqr( const f_char& SIDE
  , const f_char& TRANS, const f_int& M, const f_int& N
  , const f_int& K, const f_dbl_prec* A, const f_int& LDA
  , const f_dbl_prec* TAU, f_dbl_prec* C, const f_int& LDC
  , f_dbl_prec* WORK, const f_int& LWORK, f_int* INFO );

// DSYTRF

void dsytrf( const f_char& UPLO
  , const f_int& N, f_dbl_prec A[], const f_int& LDA
  , f_int IPIV[], f_dbl_prec WORK[], const f_int& LWORK
  , f_int* INFO );

// DSYTRS

void dsytrs( const f_char& UPLO
  , const f_int& N, const f_int& NRHS, const f_dbl_prec A[]
  , const f_int& LDA, const f_int IPIV[], f_dbl_prec B[]
  , const f_int& LDB, f_int* INFO );

// DGETRF

void dgetrf(
  const f_int& M, const f_int& N, f_dbl_prec A[], const f_int& LDA
  ,f_int IPIV[], f_int* INFO
  );

// DGETRS

void dgetrs(
  const f_char& TRANS
  ,const f_int& N, const f_int& NRHS, const f_dbl_prec A[]
  ,const f_int& LDA, const f_int IPIV[], f_dbl_prec B[]
  ,const f_int& LDB, f_int* INFO
  );


}	// end namespace LAPACK_C_Decl

#endif	// LAPACK_C_DECL_H
