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
#include <complex>


extern "C" void F77_BLAS_MANGLE(clarnv, CLARNV)
  (const int* const IDIST,
   int ISEED[],
   const int* const N,
   std::complex<float> X[]);

extern "C" void F77_BLAS_MANGLE(cpotri, CPOTRI)
  (const char* const UPLO,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(cpotrf, CPOTRF)
  (const char* const UPLO,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(cpotrs, CPOTRS)
  (const char* const UPLO,
   const int* const N,
   const int* const NRHS,
   const std::complex<float> A[],
   const int* const LDA,
   std::complex<float> B[],
   const int* const LDB,
   int* const INFO);

#ifdef HAVE_LAPACK_CLARFGP
extern "C" void F77_BLAS_MANGLE(clarfgp,CLARFGP)
  (const int* const N,                 // IN
   std::complex<float>* const ALPHA,   // IN/OUT
   std::complex<float> X[],            // IN/OUT
   const int* const INCX,              // IN
   std::complex<float>* const TAU);    // OUT
#else
#  ifdef HAVE_LAPACK_CLARFP
extern "C" void F77_BLAS_MANGLE(clarfp,CLARFP)
  (const int* const N,                 // IN
   std::complex<float>* const ALPHA,   // IN/OUT
   std::complex<float> X[],            // IN/OUT
   const int* const INCX,              // IN
   std::complex<float>* const TAU);    // OUT
#  else
extern "C" void F77_BLAS_MANGLE(clarfg,CLARFG)
  (const int* const N,                 // IN
   std::complex<float>* const ALPHA,   // IN/OUT
   std::complex<float> X[],            // IN/OUT
   const int* const INCX,              // IN
   std::complex<float>* const TAU);    // OUT
#  endif // HAVE_LAPACK_CLARFP
#endif // HAVE_LAPACK_CLARFGP

extern "C" void F77_BLAS_MANGLE(cgeqrf, CGEQRF)
  (const int* const M,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   const int* const LWORK,
   int* const INFO);

#ifdef HAVE_LAPACK_CGEQRFP
extern "C" void F77_BLAS_MANGLE(cgeqrfp, CGEQRFP)
  (const int* const M,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   const int* const LWORK,
   int* const INFO);
#endif // HAVE_LAPACK_CGEQRFP

extern "C" void F77_BLAS_MANGLE(cgeqr2, CGEQR2)
  (const int* const M,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   int* const INFO);

#ifdef HAVE_LAPACK_CGEQR2P
extern "C" void F77_BLAS_MANGLE(cgeqr2p, CGEQR2P)
  (const int* const M,
   const int* const N,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   int* const INFO);
#endif // HAVE_LAPACK_CGEQR2P

// In the complex case, Q is called UNitary rather than ORthogonal.
// This is why we have ZUNGQR and CUNGQR, rather than ZORGQR and
// CORGQR.  The interface is exactly the same as in the real case,
// though, so our LAPACK::ORMQR(), etc. wrappers have the same name
// for both the real and the complex cases.

extern "C" void F77_BLAS_MANGLE(cungqr, CUNGQR)
  (const int* const M,
   const int* const N,
   const int* const K,
   std::complex<float> A[],
   const int* const LDA,
   std::complex<float> TAU[],
   std::complex<float> WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(cunmqr, CUNMQR)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const std::complex<float> A[],
   const int* const LDA,
   const std::complex<float> TAU[],
   std::complex<float> C[],
   const int* const LDC,
   std::complex<float> WORK[],
   const int* const LWORK,
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(cunm2r, CUNM2R)
  (const char* const SIDE,
   const char* const TRANS,
   const int* const M,
   const int* const N,
   const int* const K,
   const std::complex<float> A[],
   const int* const LDA,
   const std::complex<float> TAU[],
   std::complex<float> C[],
   const int* const LDC,
   std::complex<float> WORK[],
   int* const INFO);

extern "C" void F77_BLAS_MANGLE(cgesvd, CGESVD) 
  (const char* const JOBU, 
   const char* const JOBVT, 
   const int* const M, 
   const int* const N, 
   std::complex<float> A[], 
   const int* const LDA,
   float S[], 
   std::complex<float> U[], 
   const int* const LDU, 
   std::complex<float> VT[], 
   const int* const LDVT, 
   std::complex<float> work[],
   const int* const LWORK,
   float RWORK[],
   int* const INFO);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  // If _GEQRFP is available, LAPACK::GEQRF() calls it.  If _LARFP is
  // available, LAPACK::GEQRF() calls _GEQRF, which uses _LARFP.
#ifdef HAVE_LAPACK_CGEQRFP
  template <>
  bool LAPACK<int, std::complex<float> >::QR_produces_R_factor_with_nonnegative_diagonal() { return true; }
#else
#  ifdef HAVE_LAPACK_CLARFP
  template <>
  bool LAPACK<int, std::complex<float> >::QR_produces_R_factor_with_nonnegative_diagonal() { return true; }
#  else
  template <>
  bool LAPACK<int, std::complex<float> >::QR_produces_R_factor_with_nonnegative_diagonal() { return false; }
#  endif
#endif

  ////////////////////////////////////////////////////////////////////////////
  // LARFP (implemented with _LARFGP if available, else with _LARFP if
  // available, else fall back to _LARFG)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void 
  LAPACK<int, std::complex<float> >::LARFP (const int n, 
					    std::complex<float>& alpha, 
					    std::complex<float> x[], 
					    const int incx, 
					    std::complex<float>& tau)
  {
#ifdef HAVE_LAPACK_CLARFGP
    F77_BLAS_MANGLE(clarfgp,CLARFGP) (&n, &alpha, x, &incx, &tau);
#else // Don't HAVE_LAPACK_CLARFGP
#  ifdef HAVE_LAPACK_CLARFP
    F77_BLAS_MANGLE(clarfp,CLARFP) (&n, &alpha, x, &incx, &tau);
#  else
    F77_BLAS_MANGLE(clarfg,CLARFG) (&n, &alpha, x, &incx, &tau);
#  endif // HAVE_LAPACK_CLARFP
#endif // HAVE_LAPACK_CLARFGP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQRF (implemented with _GEQRFP if available, else fall back to _GEQRF)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, std::complex<float> >::GEQRF (const int m,
					    const int n, 
					    std::complex<float> A[],
					    const int lda, 
					    std::complex<float> tau[],
					    std::complex<float> work[],
					    const int lwork,
					    int* const INFO)
  {
#ifdef HAVE_LAPACK_CGEQRFP
    F77_BLAS_MANGLE(cgeqrfp, CGEQRFP) 
      (&m, &n, A, &lda, tau, work, &lwork, INFO);
#else
    F77_BLAS_MANGLE(cgeqrf, CGEQRF) 
      (&m, &n, A, &lda, tau, work, &lwork, INFO);
#endif // HAVE_LAPACK_CGEQRFP
  }

  ////////////////////////////////////////////////////////////////////////////
  // GEQR2 (implemented with _GEQR2P if available, else fall back to _GEQR2)
  ////////////////////////////////////////////////////////////////////////////
  template <>
  void
  LAPACK<int, std::complex<float> >::GEQR2 (const int m,
					    const int n, 
					    std::complex<float> A[],
					    const int lda, 
					    std::complex<float> tau[],
					    std::complex<float> work[],
					    int* const INFO)
  {
#ifdef HAVE_LAPACK_CGEQR2P
    F77_BLAS_MANGLE(cgeqr2p, CGEQR2P) (&m, &n, A, &lda, tau, work, INFO);
#else
    F77_BLAS_MANGLE(cgeqr2, CGEQR2) (&m, &n, A, &lda, tau, work, INFO);
#endif // HAVE_LAPACK_CGEQR2P
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::
  ORMQR (const char* const side,
	 const char* const trans,
	 const int m,
	 const int n,
	 const int k,
	 const std::complex<float> A[],
	 const int lda,
	 const std::complex<float> tau[],
	 std::complex<float> C[],
	 const int ldc,
	 std::complex<float> work[],
	 const int lwork,
	 int* const INFO)
  {
    F77_BLAS_MANGLE(cunmqr, CUNMQR) 
      (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::
  ORM2R (const char* const side,
	 const char* const trans,
	 const int m,
	 const int n,
	 const int k,
	 const std::complex<float> A[],
	 const int lda,
	 const std::complex<float> tau[],
	 std::complex<float> C[],
	 const int ldc,
	 std::complex<float> work[],
	 int* const INFO)
  {
    F77_BLAS_MANGLE(cunm2r, CUNM2R) 
      (side, trans, &m, &n, &k, A, &lda, tau, C, &ldc, work, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::
  ORGQR (const int m,
	 const int n,
	 const int k,
	 std::complex<float> A[],
	 const int lda,
	 std::complex<float> tau[],
	 std::complex<float> work[],
	 const int lwork,
	 int* const INFO)
  {
    F77_BLAS_MANGLE(cungqr, CUNGQR) 
      (&m, &n, &k, A, &lda, tau, work, &lwork, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::POTRF (const char* const uplo,
					    const int n,
					    std::complex<float> A[],
					    const int lda,
					    int* const INFO)
  {
    F77_BLAS_MANGLE(cpotrf, CPOTRF) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::POTRS (const char* const uplo,
					    const int n,
					    const int nrhs,
					    const std::complex<float> A[],
					    const int lda,
					    std::complex<float> B[],
					    const int ldb,
					    int* const INFO)
  {
    F77_BLAS_MANGLE(cpotrs, CPOTRS) (uplo, &n, &nrhs, A, &lda, B, &ldb, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::POTRI (const char* const uplo, 
					    const int n, 
					    std::complex<float> A[], 
					    const int lda, 
					    int* const INFO)
  {
    F77_BLAS_MANGLE(cpotri, CPOTRI) (uplo, &n, A, &lda, INFO);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::LARNV (const int idist, 
					    int iseed[],
					    const int n,
					    std::complex<float> x[])
  {
    F77_BLAS_MANGLE(clarnv, CLARNV) (&idist, iseed, &n, x);
  }

  template <>
  void
  LAPACK<int, std::complex<float> >::
  GESVD (const char* const jobu,
	 const char* const jobvt,
	 const int m,
	 const int n,
	 std::complex<float> A[],
	 const int lda,
	 float s[],
	 std::complex<float> U[],
	 const int ldu,
	 std::complex<float> VT[],
	 const int ldvt,
	 std::complex<float> work[],
	 const int lwork,
	 float rwork[],
	 int* const INFO)
  {
    F77_BLAS_MANGLE(cgesvd, CGESVD) (jobu, jobvt, &m, &n, 
				     A, &lda, s, 
				     U, &ldu, VT, &ldvt, 
				     work, &lwork, rwork, INFO);
  }


} // namespace TSQR
