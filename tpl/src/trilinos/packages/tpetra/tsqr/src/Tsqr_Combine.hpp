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

/// \file Tsqr_Combine.hpp
/// \brief TSQR's six computational kernels.
///
#ifndef __TSQR_Combine_hpp
#define __TSQR_Combine_hpp

#include <Teuchos_ScalarTraits.hpp>
#include <Tsqr_ApplyType.hpp>
#include <Tsqr_CombineNative.hpp>


namespace TSQR {

  /// \class Combine
  /// \brief TSQR's six computational kernels
  /// \author Mark Hoemmen
  ///
  /// This class provides the six computational primitives required by
  /// TSQR.  The primitives are as follows, in which R, R_1, and R_2
  /// each represent an n x n upper triangular matrix, A represents an
  /// m x n cache block, and C_1 and C_2 represent cache blocks with
  /// some number of columns p:
  /// - Factor A (factor_first)
  /// - Apply Q factor of A to C (apply_first)
  /// - Factor [R; A] (factor_inner)
  /// - Factor [R_1; R_2] (factor_pair)
  /// - Apply Q factor of [R; A] to [C_1; C_2] (apply_inner)
  /// - Apply Q factor of [R_1; R_2] to [C_1; C_2] (apply_pair)
  ///
  /// \tparam Ordinal Type of indices into matrices.
  /// \tparam Scalar Type of entries of matrices.
  /// \tparam CombineImpl Type of a particular implementation of
  ///   Combine.  Its public interface must contain this class'
  ///   interface.
  /// 
  /// All Combine methods are implemented using CombineImpl methods
  /// with the same name.  TSQR includes three implementations of the
  /// CombineImpl interface:
  /// - \c CombineDefault, which uses LAPACK and copies in and out of
  ///   scratch space that it owns, 
  /// - \c CombineNative, a C++ in-place (no scratch space) generic 
  ///   implementation), and 
  /// - \c CombineFortran, a Fortran 9x in-place implementation for
  ///   LAPACK's four data types S, D, C, Z.
  ///
  /// The default CombineImpl is \c CombineNative, since that should
  /// work for any Ordinal and Scalar types for which LAPACK<Ordinal,
  /// Scalar> and BLAS<Ordinal, Scalar> are implemented.
  ///
  template< class Ordinal, 
	    class Scalar, 
	    class CombineImpl = CombineNative<Ordinal, Scalar, Teuchos::ScalarTraits<Scalar >::isComplex> >
  class Combine {
  public:
    /// \typedef scalar_type
    /// \brief Type of matrix entries.
    typedef Scalar scalar_type;
    /// \typedef ordinal_type
    /// \brief Type of (intranode) matrix indices.
    typedef Ordinal ordinal_type;
    /// \typedef combine_impl_type
    /// \brief Type of the implementation of Combine.
    typedef CombineImpl combine_impl_type;

    //! Constructor.
    Combine () {}

    /// Whether or not the QR factorizations computed by methods of
    /// this class produce an R factor with all nonnegative diagonal
    /// entries.  
    static bool QR_produces_R_factor_with_nonnegative_diagonal() {
      return combine_impl_type::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    /// \brief Factor the first cache block.
    ///
    /// Compute the QR factorization of the nrows by ncols matrix A
    /// (with leading dimension lda).  Overwrite the upper triangle of
    /// A with the resulting R factor, and the lower trapezoid of A
    /// (along with the length ncols tau array) with the implicitly
    /// stored Q factor.  
    ///
    /// \param nrows [in] Number of rows in A
    /// \param ncols [in] Number of columns in A
    /// \param A [in/out] On input: the nrows by ncols matrix (in
    ///   column-major order, with leading dimension lda) to factor.
    ///   On output: upper triangle contains the R factor, and lower
    ///   part contains the implicitly stored Q factor.
    /// \param lda [in] Leading dimension of A
    /// \param tau [out] Array of length ncols; on output, the 
    ///   scaling factors for the Householder reflectors
    /// \param work [out] Workspace array of length ncols
    void
    factor_first (const Ordinal nrows,
		  const Ordinal ncols,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      return impl_.factor_first (nrows, ncols, A, lda, tau, work);
    }

    /// \brief Apply the result of \c factor_first().
    ///
    /// Apply the Q factor, as computed by factor_first() and stored
    /// implicitly in A and tau, to the matrix C.
    void
    apply_first (const ApplyType& applyType,
		 const Ordinal nrows,
		 const Ordinal ncols_C,
		 const Ordinal ncols_A,
		 const Scalar A[],
		 const Ordinal lda,
		 const Scalar tau[],
		 Scalar C[],
		 const Ordinal ldc,
		 Scalar work[]) const
    {
      return impl_.apply_first (applyType, nrows, ncols_C, ncols_A, 
				A, lda, tau, C, ldc, work);
    }

    /// Apply the result of \c factor_inner().
    ///
    /// Apply the Q factor stored in [R; A] to [C_top; C_bot].  The C
    /// blocks are allowed, but not required, to have different leading
    /// dimensions (ldc_top resp. ldc_bottom).  R is upper triangular, so
    /// we do not need it; the Householder reflectors representing the Q
    /// factor are stored compactly in A (specifically, in all of A, not
    /// just the lower triangle).
    ///
    /// In the "sequential under parallel" version of TSQR, this function
    /// belongs to the sequential part (i.e., operating on cache blocks on
    /// a single processor).
    ///
    /// \param apply_type [in] NoTranspose means apply Q, Transpose
    ///   means apply Q^T, and ConjugateTranspose means apply Q^H.
    /// \param m [in]         number of rows of A
    /// \param ncols_C [in]   number of columns of [C_top; C_bot]
    /// \param ncols_Q [in]   number of columns of [R; A]
    /// \param A [in] m by ncols_Q matrix, in which the Householder 
    ///   reflectors representing the Q factor are stored
    /// \param lda [in]       leading dimension of A
    /// \param tau [in] array of length ncols_Q, storing the scaling
    ///   factors for the Householder reflectors representing Q 
    /// \param C_top [inout]  ncols_Q by ncols_C matrix
    /// \param ldc_top [in]   leading dimension of C_top
    /// \param C_bot [inout]  m by ncols_C matrix
    /// \param ldc_bot [in]   leading dimension of C_bot
    /// \param work [out]     workspace array of length ncols_C
    void
    apply_inner (const ApplyType& apply_type,
		 const Ordinal m,
		 const Ordinal ncols_C,
		 const Ordinal ncols_Q,
		 const Scalar A[],
		 const Ordinal lda,
		 const Scalar tau[],
		 Scalar C_top[],
		 const Ordinal ldc_top,
		 Scalar C_bot[],
		 const Ordinal ldc_bot,
		 Scalar work[]) const
    {
      impl_.apply_inner (apply_type, m, ncols_C, ncols_Q, 
			 A, lda, tau, 
			 C_top, ldc_top, C_bot, ldc_bot, work);
    }

    /// \brief Factor [R; A] for square upper triangular R and cache block A.
    ///
    /// Perform one "inner" QR factorization step of sequential / parallel
    /// TSQR.  (In either case, only one processor calls this routine.)
    ///
    /// In the "sequential under parallel" version of TSQR, this function
    /// belongs to the sequential part (i.e., operating on cache blocks on
    /// a single processor).  Only the first cache block $A_0$ is factored
    /// as $Q_0 R_0 = A_0$ (see tsqr_factor_first); subsequent cache blocks
    /// $A_k$ are factored using this routine, which combines them with 
    /// $R_{k-1}$.
    ///
    /// Here is the matrix to factor:
    /// \[
    /// \begin{pmatrix}
    /// R_{k-1} \\      % $A_{k-1}$ is $m_{k-1} \times n$ with $m_{k-1} \geq n$
    /// A_k     \\      % $m_k \times n$ with $m_k \geq n$
    /// \end{pmatrix}
    /// \]
    ///
    /// Since $R_{k-1}$ is n by n upper triangular, we can store the
    /// Householder reflectors (representing the Q factor of $[R_{k-1};
    /// A_k]$) entirely in $A_k$ (specifically, in all of $A_k$, not just
    /// below the diagonal).
    ///
    /// \param m [in] Number of rows in the "bottom" block to factor.  
    ///   The number of rows in the top block doesn't matter, given the 
    ///   assumptions above, as long as $m_{k-1} \geq n$.
    /// \param n [in] Number of columns (same in both blocks)
    /// \param R [inout] "Top" upper triangular n by n block $R_{k-1}$.
    ///   Overwritten with the new R factor $R_k$ of $[R_{k-1}; A_k]$.
    /// \param ldr [in] Leading dimension of R
    /// \param A [inout] "Bottom" dense m by n block $A_k$.  Overwritten 
    ///   with the Householder reflectors representing the Q factor of 
    ///   $[R_{k-1}; A_k]$.
    /// \param tau [out] Scaling factors of the Householder reflectors.  
    ///   Corresponds to the TAU output of LAPACK's _GEQRF.
    /// \param work [out] Workspace (length >= n; don't need lwork or 
    ///   workspace query)
    void
    factor_inner (const Ordinal m,
		  const Ordinal n,
		  Scalar R[],
		  const Ordinal ldr,
		  Scalar A[],
		  const Ordinal lda,
		  Scalar tau[],
		  Scalar work[]) const
    {
      impl_.factor_inner (m, n, R, ldr, A, lda, tau, work);
    }

    /// \brief Factor the pair of square upper triangular matrices [R_top; R_bot].  
    ///
    /// Store the resulting R factor in R_top, and the resulting
    /// Householder reflectors implicitly in R_bot and tau.
    ///
    /// \param n [in]         Number of rows and columns of each of R_top and R_bot
    /// \param R_top [inout]  n by n upper triangular matrix 
    /// \param ldr_top [in]   Leading dimension of R_top
    /// \param R_bot [inout]  n by n upper triangular matrix 
    /// \param ldr_bot [in]   Leading dimension of R_bot
    /// \param tau [out]      Scaling factors for Householder reflectors
    /// \param work [out]     Workspace array (of length >= n)
    ///
    void
    factor_pair (const Ordinal n,
		 Scalar R_top[],
		 const Ordinal ldr_top,
		 Scalar R_bot[],
		 const Ordinal ldr_bot,
		 Scalar tau[],
		 Scalar work[]) const
    {
      impl_.factor_pair (n, R_top, ldr_top, R_bot, ldr_bot, tau, work);      
    }

    /// \brief Apply the result of \c factor_pair().
    ///
    /// Apply Q factor (or Q^T or Q^H) of the 2*ncols_Q by ncols_Q
    /// matrix [R_top; R_bot] (stored in R_bot and tau) to the
    /// 2*ncols_Q by ncols_C matrix [C_top; C_bot].  The two blocks
    /// C_top and C_bot may have different leading dimensions (ldc_top
    /// resp. ldc_bot).
    ///
    /// \param apply_type [in] NoTranspose means apply Q, Transpose
    ///   means apply Q^T, and ConjugateTranspose means apply Q^H.
    void
    apply_pair (const ApplyType& apply_type,
		const Ordinal ncols_C, 
		const Ordinal ncols_Q, 
		const Scalar R_bot[], 
		const Ordinal ldr_bot,
		const Scalar tau[], 
		Scalar C_top[], 
		const Ordinal ldc_top, 
		Scalar C_bot[], 
		const Ordinal ldc_bot, 
		Scalar work[]) const
    {
      impl_.apply_pair (apply_type, ncols_C, ncols_Q, 
			R_bot, ldr_bot, tau, 
			C_top, ldc_top, C_bot, ldc_bot, work);
    }

  private:
    //! The implementation of Combine.
    combine_impl_type impl_;
  };

} // namespace TSQR

#endif // __TSQR_Combine_hpp
