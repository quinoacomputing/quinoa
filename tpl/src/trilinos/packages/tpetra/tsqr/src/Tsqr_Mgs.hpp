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

#ifndef __TSQR_Tsqr_Mgs_hpp
#define __TSQR_Tsqr_Mgs_hpp

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility> // std::pair

#include <Tsqr_MessengerBase.hpp>
#include <Tsqr_Util.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>

// #define MGS_DEBUG 1
#ifdef MGS_DEBUG
#  include <iostream>
using std::cerr;
using std::endl;
#endif // MGS_DEBUG


namespace TSQR {

  /// \class MGS
  /// \brief Distributed-memory parallel implementation of Modified Gram-Schmidt.
  template<class LocalOrdinal, class Scalar>
  class MGS {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal ordinal_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;

    /// \brief Constructor
    ///
    /// \param messenger [in/out] Communicator wrapper instance.
    ///
    MGS (const Teuchos::RCP< MessengerBase< Scalar > >& messenger) :
      messenger_ (messenger) {}

    /// \brief Does the R factor have a nonnegative diagonal?
    ///
    /// MGS implements a QR factorization (of a distributed matrix).
    /// Some, but not all, QR factorizations produce an R factor whose
    /// diagonal may include negative entries.  This Boolean tells you
    /// whether MGS promises to compute an R factor whose diagonal
    /// entries are all nonnegative.
    ///
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return true;
    }

    //! Use Modified Gram-Schmidt to orthogonalize a matrix A in place.
    void
    mgs (const LocalOrdinal nrows_local,
         const LocalOrdinal ncols,
         Scalar A_local[],
         const LocalOrdinal lda_local,
         Scalar R[],
         const LocalOrdinal ldr);

  private:
    Teuchos::RCP<MessengerBase<Scalar> > messenger_;
  };


  namespace details {

    template<class LocalOrdinal, class Scalar>
    class MgsOps {
    public:
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef typename STS::magnitudeType magnitude_type;

      MgsOps (const Teuchos::RCP< MessengerBase< Scalar > >& messenger) :
        messenger_ (messenger) {}

      void
      axpy (const LocalOrdinal nrows_local,
            const Scalar alpha,
            const Scalar x_local[],
            Scalar y_local[]) const
      {
        for (LocalOrdinal i = 0; i < nrows_local; ++i)
          y_local[i] = y_local[i] + alpha * x_local[i];
      }

      void
      scale (const LocalOrdinal nrows_local,
             Scalar x_local[],
             const Scalar denom) const
      {
        for (LocalOrdinal i = 0; i < nrows_local; ++i)
          x_local[i] = x_local[i] / denom;
      }

      /// $y^* \cdot x$: conjugate transpose when Scalar is complex,
      /// else regular transpose.
      Scalar
      dot (const LocalOrdinal nrows_local,
           const Scalar x_local[],
           const Scalar y_local[])
      {
        Scalar local_result (0);

#ifdef MGS_DEBUG
        // for (LocalOrdinal k = 0; k != nrows_local; ++k)
        //   cerr << "(x[" << k << "], y[" << k << "]) = (" << x_local[k] << "," << y_local[k] << ")" << " ";
        //   cerr << endl;
#endif // MGS_DEBUG

        for (LocalOrdinal i = 0; i < nrows_local; ++i)
          local_result += x_local[i] * STS::conjugate (y_local[i]);

#ifdef MGS_DEBUG
          // cerr << "-- Final value on this proc = " << local_result << endl;
#endif // MGS_DEBUG

        // FIXME (mfh 23 Apr 2010) Does MPI_SUM do the right thing for
        // complex or otherwise general MPI data types?  Perhaps an MPI_Op
        // should belong in the MessengerBase...
        return messenger_->globalSum (local_result);
      }

      magnitude_type
      norm2 (const LocalOrdinal nrows_local,
             const Scalar x_local[])
      {
        Scalar localResult (0);

        // Doing the right thing in the complex case requires taking
        // an absolute value.  We want to avoid this additional cost
        // in the real case, which is why we check is_complex.
        if (STS::isComplex)
          {
            for (LocalOrdinal i = 0; i < nrows_local; ++i)
              {
                const Scalar xi = STS::magnitude (x_local[i]);
                localResult += xi * xi;
              }
          }
        else
          {
            for (LocalOrdinal i = 0; i < nrows_local; ++i)
              {
                const Scalar xi = x_local[i];
                localResult += xi * xi;
              }
          }
        const Scalar globalResult = messenger_->globalSum (localResult);
        // sqrt doesn't make sense if the type of Scalar is complex,
        // even if the imaginary part of global_result is zero.
        return STS::squareroot (STS::magnitude (globalResult));
      }

      Scalar
      project (const LocalOrdinal nrows_local,
               const Scalar q_local[],
               Scalar v_local[])
      {
        const Scalar coeff = this->dot (nrows_local, v_local, q_local);
        this->axpy (nrows_local, -coeff, q_local, v_local);
        return coeff;
      }

    private:
      Teuchos::RCP< MessengerBase< Scalar > > messenger_;
    };
  } // namespace details


  template<class LocalOrdinal, class Scalar>
  void
  MGS<LocalOrdinal, Scalar>::mgs (const LocalOrdinal nrows_local,
                                  const LocalOrdinal ncols,
                                  Scalar A_local[],
                                  const LocalOrdinal lda_local,
                                  Scalar R[],
                                  const LocalOrdinal ldr)
  {
    details::MgsOps<LocalOrdinal, Scalar> ops (messenger_);

    for (LocalOrdinal j = 0; j < ncols; ++j)
      {
        Scalar* const v = &A_local[j*lda_local];
        for (LocalOrdinal i = 0; i < j; ++i)
          {
            const Scalar* const q = &A_local[i*lda_local];
            R[i + j*ldr] = ops.project (nrows_local, q, v);
#ifdef MGS_DEBUG
            if (my_rank == 0)
              cerr << "(i,j) = (" << i << "," << j << "): coeff = " << R[i + j*ldr] << endl;
#endif // MGS_DEBUG
          }
        const magnitude_type denom = ops.norm2 (nrows_local, v);
#ifdef MGS_DEBUG
          if (my_rank == 0)
            cerr << "j = " << j << ": denom = " << denom << endl;
#endif // MGS_DEBUG

        // FIXME (mfh 29 Apr 2010)
        //
        // NOTE IMPLICIT CAST.  This should work for complex numbers.
        // If it doesn't work for your Scalar data type, it means that
        // you need a different data type for the diagonal elements of
        // the R factor, than you need for the other elements.  This
        // is unlikely if we're comparing MGS against a Householder QR
        // factorization; I don't really understand how the latter
        // would work (not that it couldn't be given a sensible
        // interpretation) in the case of Scalars that aren't plain
        // old real or complex numbers.
        R[j + j*ldr] = Scalar (denom);
        ops.scale (nrows_local, v, denom);
      }
  }

} // namespace TSQR

#endif // __TSQR_Tsqr_Mgs_hpp
