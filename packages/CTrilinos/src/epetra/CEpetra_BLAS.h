#ifndef CEPETRA_BLAS_H
#define CEPETRA_BLAS_H

/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


#include "CTrilinos_config.h"


/*! @file CEpetra_BLAS.h
 * @brief Wrappers for Epetra_BLAS */

/* True C header file! */

#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif



/*! @name ID struct conversion functions */
/*@{*/

/*! @brief Changes the ID struct from the universal
   (generalized) struct type to the class-specific one.
*/
CT_Epetra_BLAS_ID_t Epetra_BLAS_Degeneralize ( 
  CTrilinos_Universal_ID_t id );

/*! @brief Changes the ID struct from the class-specific
   struct type to the universal (generalized) one.
*/
CTrilinos_Universal_ID_t Epetra_BLAS_Generalize ( 
  CT_Epetra_BLAS_ID_t id );

/*@}*/

/*! @name Epetra_BLAS constructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   Epetra_BLAS::Epetra_BLAS(void)
*/
CT_Epetra_BLAS_ID_t Epetra_BLAS_Create (  );

/*! @brief Wrapper for 
   Epetra_BLAS::Epetra_BLAS(const Epetra_BLAS& BLAS)
*/
CT_Epetra_BLAS_ID_t Epetra_BLAS_Duplicate ( 
  CT_Epetra_BLAS_ID_t BLASID );

/*@}*/

/*! @name Epetra_BLAS destructor wrappers */
/*@{*/

/*! @brief Wrapper for 
   virtual Epetra_BLAS::~Epetra_BLAS(void)
*/
void Epetra_BLAS_Destroy ( CT_Epetra_BLAS_ID_t * selfID );

/*@}*/

/*! @name Epetra_BLAS member wrappers */
/*@{*/

/*! @brief Wrapper for 
   float Epetra_BLAS::ASUM(const int N, const float * X, const int INCX = 1) const
*/
float Epetra_BLAS_ASUM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX );

/*! @brief Wrapper for 
   double Epetra_BLAS::ASUM(const int N, const double * X, const int INCX = 1) const
*/
double Epetra_BLAS_ASUM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX );

/*! @brief Wrapper for 
   float Epetra_BLAS::DOT(const int N, const float * X, const float * Y, const int INCX = 1, const int INCY = 1) const
*/
float Epetra_BLAS_DOT_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const float * Y, const int INCX, const int INCY );

/*! @brief Wrapper for 
   double Epetra_BLAS::DOT(const int N, const double * X, const double * Y, const int INCX = 1, const int INCY = 1) const
*/
double Epetra_BLAS_DOT_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const double * Y, const int INCX, const int INCY );

/*! @brief Wrapper for 
   float Epetra_BLAS::NRM2(const int N, const float * X, const int INCX = 1) const
*/
float Epetra_BLAS_NRM2_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX );

/*! @brief Wrapper for 
   double Epetra_BLAS::NRM2(const int N, const double * X, const int INCX = 1) const
*/
double Epetra_BLAS_NRM2_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX );

/*! @brief Wrapper for 
   void Epetra_BLAS::SCAL( const int N, const float ALPHA, float * X, const int INCX = 1) const
*/
void Epetra_BLAS_SCAL_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  float * X, const int INCX );

/*! @brief Wrapper for 
   void Epetra_BLAS::SCAL( const int N, const double ALPHA, double * X, const int INCX = 1) const
*/
void Epetra_BLAS_SCAL_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  double * X, const int INCX );

/*! @brief Wrapper for 
   void Epetra_BLAS::COPY( const int N, const float * X, float * Y, const int INCX = 1, const int INCY = 1) const
*/
void Epetra_BLAS_COPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  float * Y, const int INCX, const int INCY );

/*! @brief Wrapper for 
   void Epetra_BLAS::COPY( const int N, const double * X, double * Y, const int INCX = 1, const int INCY = 1) const
*/
void Epetra_BLAS_COPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  double * Y, const int INCX, const int INCY );

/*! @brief Wrapper for 
   int Epetra_BLAS::IAMAX( const int N, const float * X, const int INCX = 1) const
*/
int Epetra_BLAS_IAMAX_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX );

/*! @brief Wrapper for 
   int Epetra_BLAS::IAMAX( const int N, const double * X, const int INCX = 1) const
*/
int Epetra_BLAS_IAMAX_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX );

/*! @brief Wrapper for 
   void Epetra_BLAS::AXPY( const int N, const float ALPHA, const float * X, float * Y, const int INCX = 1, const int INCY = 1) const
*/
void Epetra_BLAS_AXPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  const float * X, float * Y, const int INCX, const int INCY );

/*! @brief Wrapper for 
   void Epetra_BLAS::AXPY( const int N, const double ALPHA, const double * X, double * Y, const int INCX = 1, const int INCY = 1) const
*/
void Epetra_BLAS_AXPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  const double * X, double * Y, const int INCX, const int INCY );

/*! @brief Wrapper for 
   void Epetra_BLAS::GEMV(const char TRANS, const int M, const int N, const float ALPHA, const float * A, const int LDA, const float * X, const float BETA, float * Y, const int INCX = 1, const int INCY = 1) const
*/
void Epetra_BLAS_GEMV_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const float ALPHA, const float * A, const int LDA, 
  const float * X, const float BETA, float * Y, const int INCX, 
  const int INCY );

/*! @brief Wrapper for 
   void Epetra_BLAS::GEMV(const char TRANS, const int M, const int N, const double ALPHA, const double * A, const int LDA, const double * X, const double BETA, double * Y, const int INCX = 1, const int INCY = 1) const
*/
void Epetra_BLAS_GEMV_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const double ALPHA, const double * A, const int LDA, 
  const double * X, const double BETA, double * Y, const int INCX, 
  const int INCY );

/*! @brief Wrapper for 
   void Epetra_BLAS::GEMM(const char TRANSA, const char TRANSB, const int M, const int N, const int K, const float ALPHA, const float * A, const int LDA, const float * B, const int LDB, const float BETA, float * C, const int LDC) const
*/
void Epetra_BLAS_GEMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const float ALPHA, 
  const float * A, const int LDA, const float * B, const int LDB, 
  const float BETA, float * C, const int LDC );

/*! @brief Wrapper for 
   void Epetra_BLAS::GEMM(const char TRANSA, const char TRANSB, const int M, const int N, const int K, const double ALPHA, const double * A, const int LDA, const double * B, const int LDB, const double BETA, double * C, const int LDC) const
*/
void Epetra_BLAS_GEMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const double ALPHA, 
  const double * A, const int LDA, const double * B, const int LDB, 
  const double BETA, double * C, const int LDC );

/*! @brief Wrapper for 
   void Epetra_BLAS::SYMM(const char SIDE, const char UPLO, const int M, const int N, const float ALPHA, const float * A, const int LDA, const float * B, const int LDB, const float BETA, float * C, const int LDC) const
*/
void Epetra_BLAS_SYMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const float ALPHA, const float * A, 
  const int LDA, const float * B, const int LDB, const float BETA, 
  float * C, const int LDC );

/*! @brief Wrapper for 
   void Epetra_BLAS::SYMM(const char SIDE, const char UPLO, const int M, const int N, const double ALPHA, const double * A, const int LDA, const double * B, const int LDB, const double BETA, double * C, const int LDC) const
*/
void Epetra_BLAS_SYMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const double ALPHA, const double * A, 
  const int LDA, const double * B, const int LDB, const double BETA, 
  double * C, const int LDC );

/*! @brief Wrapper for 
   void Epetra_BLAS::TRMM(const char SIDE, const char UPLO, const char TRANSA, const char DIAG, const int M, const int N, const float ALPHA, const float * A, const int LDA, float * B, const int LDB) const
*/
void Epetra_BLAS_TRMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const float ALPHA, const float * A, const int LDA, float * B, 
  const int LDB );

/*! @brief Wrapper for 
   void Epetra_BLAS::TRMM(const char SIDE, const char UPLO, const char TRANSA, const char DIAG, const int M, const int N, const double ALPHA, const double * A, const int LDA, double * B, const int LDB) const
*/
void Epetra_BLAS_TRMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const double ALPHA, const double * A, const int LDA, double * B, 
  const int LDB );

/*@}*/


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_BLAS_H */

