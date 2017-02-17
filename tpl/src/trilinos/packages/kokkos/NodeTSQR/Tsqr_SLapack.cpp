//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Tsqr_Lapack.hpp>


extern "C" void F77_BLAS_MANGLE(slarnv, SLARNV)
  (const int* const IDIST,
   int ISEED[],
   const int* const N,
   float X[]);

extern "C" void F77_BLAS_MANGLE(spotri, SPOTRI)
  (const char* const UPLO,
   const int* const N,
   float A[],
   const int* const LDA,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(spotrf, SPOTRF)
  (const char* const UPLO,
   const int* const N,
   float A[],
   const int* const LDA,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(spotrs, SPOTRS)
  (const char* const UPLO,
   const int* const N,
   const int* const NRHS,
   const float A[],
   const int* const LDA,
   float B[],
   const int* const LDB,
   int* const INFO);

#ifdef HAVE_LAPACK_SLARFGP
extern "C" void F77_BLAS_MANGLE(slarfgp,SLARFGP)
  (const int* const N,    // IN
   float* const ALPHA,   // IN/OUT
   float X[],            // IN/OUT
   const int* const INCX, // IN
   float* const TAU);    // OUT
#else
#  ifdef HAVE_LAPACK_SLARFP
extern "C" void F77_BLAS_MANGLE(slarfp,SLARFP)
  (const int* const N,    // IN
   float* const ALPHA,   // IN/OUT
   float X[],            // IN/OUT
   const int* const INCX, // IN
   float* const TAU);    // OUT
#  else
extern "C" void F77_BLAS_MANGLE(slarfg,SLARFG)
  (const int* const N,    // IN
   float* const ALPHA,   // IN/OUT
   float X[],            // IN/OUT
   const int* const INCX, // IN
   float* const TAU);    // OUT
#  endif // HAVE_LAPACK_SLARFP
#endif // HAVE_LAPACK_SLARFGP

extern "C" void F77_BLAS_MANGLE(sgeqrf, SGEQRF)
  (const int* const M,
   const int* const N,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   const int* const LWORK,
   int* const INFO);

#ifdef HAVE_LAPACK_SGEQRFP
extern "C" void F77_BLAS_MANGLE(sgeqrfp, SGEQRFP)
  (const int* const M,
   const int* const N,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   const int* const LWORK,
   int* const INFO);
#endif // HAVE_LAPACK_SGEQRFP

extern "C" void F77_BLAS_MANGLE(sgeqr2, SGEQR2)
  (const int* const M,
   const int* const N,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   int* const INFO);

#ifdef HAVE_LAPACK_SGEQR2P
extern "C" void F77_BLAS_MANGLE(sgeqr2p, SGEQR2P)
  (const int* const M,
   const int* const N,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   int* const INFO);
#endif // HAVE_LAPACK_SGEQR2P

extern "C" void F77_BLAS_MANGLE(sormqr, SORMQR)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const float A[],
   const int* const LDA,
   const float TAU[],
   float C[],
   const int* const LDC,
   float WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(sorm2r, SORM2R)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const float A[],
   const int* const LDA,
   const float TAU[],
   float C[],
   const int* const LDC,
   float WORK[],
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(sorgqr, SORGQR)
  (const int* const M,
   const int* const N,
   const int* const K,
   float A[],
   const int* const LDA,
   float TAU[],
   float WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(sgesvd, SGESVD) 
  (const char* const JOBU, 
   const char* const JOBVT, 
   const int* const M, 
   const int* const N, 
   float A[], 
   const int* const LDA,
   float S[], 
   float U[], 
   const int* const LDU, 
   float VT[], 
   const int* const LDVT, 
   float work[],
   const int* const LWORK,
   float RWORK[],
   int* const INFO);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // If _GEQRFP is available, LAPACK::GEQRF() calls it.  If _LARFP is
  // available, LAPACK::GEQRF() calls _GEQRF, which uses _LARFP.
#ifdef HAVE_LAPACK_SGEQRFP
  template <>
  bool LAPACK<int, float >::QR_produces_R_factor_with_nonnegative_diagonal() { return true; }
#else // Don't HAVE_LAPACK_SGEQRFP
#  ifdef HAVE_LAPACK_SLARFP
  template <>
  bool LAPACK<int, float >::QR_produces_R_factor_with_nonnegative_diagonal() { return true; }
#  else
  template <>
  bool LAPACK<int, float >::QR_produces_R_factor_with_nonnegative_diagonal() { return false; }
#  endif
#endif

  ////////////////////////////////////////////////////////////////////////////
  // LARFP (implemented with _LARFGP if available, else with _LARFP if
  // available, else fall back to _LARFG)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void 
  LAPACK<int, float >::LARFP (const int n, 
			      float& alpha, 
			      float x[], 
			      const int incx, 
			      float& tau)
  {
#ifdef HAVE_LAPACK_SLARFGP
    F77_BLAS_MANGLE(slarfgp,SLARFGP) (&n, &alpha, x, &incx, &tau);
#else // Don't HAVE_LAPACK_SLARFGP
#  ifdef HAVE_LAPACK_SLARFP
    F77_BLAS_MANGLE(slarfp,SLARFP) (&n, &alpha, x, &incx, &tau);
#  else
    F77_BLAS_MANGLE(slarfg,SLARFG) (&n, &alpha, x, &incx, &tau);
#  endif // HAVE_LAPACK_SLARFP
#endif // HAVE_LAPACK_SLARFGP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQRF (implemented with _GEQRFP if available, else fall back to _GEQRF)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, float >::GEQRF (const int m,
			      const int n, 
			      float A[],
			      const int lda, 
			      float tau[],
			      float work[],
			      const int lwork,
			      int* const INFO)
  {
#ifdef HAVE_LAPACK_SGEQRFP
    F77_BLAS_MANGLE(sgeqrfp, SGEQRFP) 
      (&m, &n, A, &lda, tau, work, &lwork, INFO);
#else
    F77_BLAS_MANGLE(sgeqrf, SGEQRF) 
      (&m, &n, A, &lda, tau, work, &lwork, INFO);
#endif // HAVE_LAPACK_SGEQRFP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQR2 (implemented with _GEQR2P if available, else fall back to _GEQR2)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, float >::GEQR2 (const int m,
			      const int n, 
			      float A[],
			      const int lda, 
			      float tau[],
			      float work[],
			      int* const INFO)
  {
#ifdef HAVE_LAPACK_SGEQR2P
    F77_BLAS_MANGLE(sgeqr2p, SGEQR2P) (&m, &n, A, &lda, tau, work, INFO);
#else
    F77_BLAS_MANGLE(sgeqr2, SGEQR2) (&m, &n, A, &lda, tau, work, INFO);
#endif // HAVE_LAPACK_SGEQR2P
  }

  template <>
  void
  LAPACK<int, float >::ORMQR (const char* const side,
			      const char* const trans,
			      const int m,
			      const int n,
			      const int k,
			      const float A[],
			      const int lda,
			      const float tau[],
			      float C[],
			      const int ldc,
			      float work[],
			      const int lwork,
			      int* const INFO)
  {
    F77_BLAS_MANGLE(sormqr, SORMQR) 
      (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, float >::ORM2R (const char* const side,
			      const char* const trans,
			      const int m,
			      const int n,
			      const int k,
			      const float A[],
			      const int lda,
			      const float tau[],
			      float C[],
			      const int ldc,
			      float work[],
			      int* const INFO)
  {
    F77_BLAS_MANGLE(sorm2r, SORM2R) 
      (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, INFO);
  }

  template <>
  void
  LAPACK<int, float >::ORGQR (const int m,
			      const int n,
			      const int k,
			      float A[],
			      const int lda,
			      float tau[],
			      float work[],
			      const int lwork,
			      int* const INFO)
  {
    F77_BLAS_MANGLE(sorgqr, SORGQR) 
      (&m, &n, &k, A, &lda, tau, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, float >::POTRF (const char* const uplo,
			      const int n,
			      float A[],
			      const int lda,
			      int* const INFO)
  {
    F77_BLAS_MANGLE(spotrf, SPOTRF) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, float >::POTRS (const char* const uplo,
			      const int n,
			      const int nrhs,
			      const float A[],
			      const int lda,
			      float B[],
			      const int ldb,
			      int* const INFO)
  {
    F77_BLAS_MANGLE(spotrs, SPOTRS) (uplo, &n, &nrhs, A, &lda, B, &ldb, INFO);
  }

  template <>
  void
  LAPACK<int, float >::POTRI (const char* const uplo, 
			      const int n, 
			      float A[], 
			      const int lda, 
			      int* const INFO)
  {
    F77_BLAS_MANGLE(spotri, SPOTRI) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, float >::LARNV (const int idist, 
			      int iseed[],
			      const int n,
			      float x[])
  {
    F77_BLAS_MANGLE(slarnv, SLARNV) (&idist, iseed, &n, x);
  }

  template <>
  void
  LAPACK<int, float >::GESVD (const char* const jobu,
			      const char* const jobvt,
			      const int m,
			      const int n,
			      float A[],
			      const int lda,
			      float s[],
			      float U[],
			      const int ldu,
			      float VT[],
			      const int ldvt,
			      float work[],
			      const int lwork,
			      float rwork[],
			      int* const INFO)
  {
    F77_BLAS_MANGLE(sgesvd, SGESVD) (jobu, jobvt, &m, &n, 
				     A, &lda, s, 
				     U, &ldu, VT, &ldvt, 
				     work, &lwork, rwork, INFO);
  }

} // namespace TSQR
