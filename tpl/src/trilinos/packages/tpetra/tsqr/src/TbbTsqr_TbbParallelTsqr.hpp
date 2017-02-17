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

#ifndef __TSQR_TBB_TbbParallelTsqr_hpp
#define __TSQR_TBB_TbbParallelTsqr_hpp

#include <tbb/tbb.h>
#include <tbb/task_scheduler_init.h>

#include <TbbTsqr_FactorTask.hpp>
#include <TbbTsqr_ApplyTask.hpp>
#include <TbbTsqr_ExplicitQTask.hpp>
#include <TbbTsqr_RevealRankTask.hpp>
#include <TbbTsqr_CacheBlockTask.hpp>
#include <TbbTsqr_UnCacheBlockTask.hpp>
#include <TbbTsqr_FillWithZerosTask.hpp>

#include <Tsqr_ApplyType.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <algorithm>
#include <limits>


namespace TSQR {
  namespace TBB {

    /// \class TbbParallelTsqr
    /// \brief Parallel implementation of \c TbbTsqr.
    /// \author Mark Hoemmen
    ///
    /// This class implements the functionality of \c TbbTsqr.
    /// It is not meant to be seen by users of \c TbbTsqr.
    ///
    /// The third template parameter, TimerType, allows different
    /// timer implementations.  TbbParallelTsqr times each task's
    /// invocations of \c SequentialTsqr::factor() and \c
    /// SequentialTsqr::apply().  \c TrivialTimer is a "timer" that
    /// does nothing, in case you don't want to invoke timers.
    template<class LocalOrdinal, class Scalar, class TimerType>
    class TbbParallelTsqr {
    private:
      typedef MatView<LocalOrdinal, Scalar> mat_view_type;
      typedef ConstMatView<LocalOrdinal, Scalar> const_mat_view_type;
      typedef std::pair<mat_view_type, mat_view_type> split_t;
      typedef std::pair<const_mat_view_type, const_mat_view_type> const_split_t;
      typedef std::pair<const_mat_view_type, mat_view_type> top_blocks_t;
      typedef std::vector<top_blocks_t> array_top_blocks_t;

      template<class MatrixViewType>
      MatrixViewType
      top_block_helper (const size_t P_first,
                        const size_t P_last,
                        const MatrixViewType& C,
                        const bool contiguous_cache_blocks) const
      {
        if (P_first > P_last)
          throw std::logic_error ("P_first > P_last");
        else if (P_first == P_last)
          return seq_.top_block (C, contiguous_cache_blocks);
        else
          {
            typedef std::pair<MatrixViewType, MatrixViewType> split_type;

            // Divide [P_first, P_last] into two intervals: [P_first,
            // P_mid] and [P_mid+1, P_last].  Recurse on the first
            // interval [P_first, P_mid].
            const size_t P_mid = (P_first + P_last) / 2;
            split_type C_split = partitioner_.split (C, P_first, P_mid, P_last,
                                                     contiguous_cache_blocks);
            // The partitioner may decide that the current block C has
            // too few rows to be worth splitting.  In that case,
            // C_split.first should be the same block as C, and
            // C_split.second (the bottom block) will be empty.  We
            // deal with this in the same way as the base case
            // (P_first == P_last) above.
            if (C_split.second.empty() || C_split.second.nrows() == 0)
              return seq_.top_block (C_split.first, contiguous_cache_blocks);
            else
              return top_block_helper (P_first, P_mid, C_split.first,
                                       contiguous_cache_blocks);
          }
      }

    public:
      typedef Scalar scalar_type;
      typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
      typedef LocalOrdinal ordinal_type;

      /// Whether or not this QR factorization produces an R factor
      /// with all nonnegative diagonal entries.
      static bool QR_produces_R_factor_with_nonnegative_diagonal() {
        typedef Combine<LocalOrdinal, Scalar> combine_type;
        //typedef LAPACK<LocalOrdinal, Scalar> lapack_type;

        const bool combineMakesNonnegDiag =
          combine_type::QR_produces_R_factor_with_nonnegative_diagonal ();
        //const bool lapackMakesNonnegDiag =
        //  lapack_type::QR_produces_R_factor_with_nonnegative_diagonal ();
        const bool lapackMakesNonnegDiag = false;
        return combineMakesNonnegDiag && lapackMakesNonnegDiag;
      }

      /// \typedef SeqOutput
      /// \brief Results of SequentialTsqr for each core.
      typedef typename SequentialTsqr<LocalOrdinal, Scalar>::FactorOutput SeqOutput;

      /// \typedef ParOutput
      /// \brief Array of numTasks_ "local tau arrays" from parallel TSQR.
      ///
      /// (Local Q factors are stored in place.)
      typedef std::vector<std::vector<Scalar> > ParOutput;

      /// \typedef FactorOutput
      /// \brief Partial representation of the Q factor.
      ///
      /// The \c factor() method returns a pair: the results of
      /// SequentialTsqr for data on each core, and the results of
      /// combining the data on the cores.
      typedef typename std::pair<std::vector<SeqOutput>, ParOutput> FactorOutput;

      /// \brief Constructor.
      ///
      /// \param numTasks [in] Number of parallel tasks to use in the
      ///   factorization.  This should be >= the number of cores with
      ///   which Intel TBB was initialized.
      /// \param cacheSizeHint [in] Cache size hint in bytes.  Zero
      ///   means that TSQR will pick a reasonable nonzero default.
      TbbParallelTsqr (const size_t numTasks = 1,
                       const size_t cacheSizeHint = 0) :
        seq_ (cacheSizeHint),
        min_seq_factor_timing_ (std::numeric_limits<double>::max()),
        max_seq_factor_timing_ (std::numeric_limits<double>::min()),
        min_seq_apply_timing_ (std::numeric_limits<double>::max()),
        max_seq_apply_timing_ (std::numeric_limits<double>::min())
      {
        if (numTasks < 1)
          numTasks_ = 1; // default is no parallelism
        else
          numTasks_ = numTasks;
      }

      /// \brief Constructor (that takes a parameter list).
      ///
      /// \param plist [in/out] On input: list of parameters.  On
      ///   output: missing parameters are filled in with default
      ///   values.
      ///
      /// For a list of accepted parameters and thei documentation,
      /// see the parameter list returned by \c getValidParameters().
      TbbParallelTsqr (const Teuchos::RCP<Teuchos::ParameterList>& plist) :
        seq_ (plist), // SequentialTsqr has a plist-accepting constructor.
        numTasks_ (1),  // Set a safe default for now.
        min_seq_factor_timing_ (std::numeric_limits<double>::max()),
        max_seq_factor_timing_ (std::numeric_limits<double>::min()),
        min_seq_apply_timing_ (std::numeric_limits<double>::max()),
        max_seq_apply_timing_ (std::numeric_limits<double>::min())
      {
        if (! plist.is_null()) {
          const int defaultNumTasks = 1; // A reasonable safe default value.
          int numTasks = plist->get ("Num Tasks", defaultNumTasks);
          if (numTasks < 1) { // Default is no parallelism.
            plist->set ("Num Tasks", defaultNumTasks);
          }
          numTasks_ = numTasks;
        }
      }

      Teuchos::RCP<const Teuchos::ParameterList>
      getValidParameters () const
      {
        using Teuchos::ParameterList;
        using Teuchos::parameterList;
        using Teuchos::RCP;

        // TbbTsqr recursively divides the tall skinny matrix on the
        // node into TBB tasks.  Each task works on a block row.  The
        // TBB task scheduler ensures that oversubscribing TBB tasks
        // won't oversubscribe cores, so it's OK if
        // default_num_threads() is too many.  For example, TBB might
        // say default_num_threads() is the number of cores on the
        // node, but the TBB task scheduler might have been
        // initialized with the number of cores per NUMA region, for
        // hybrid MPI + TBB parallelism.
        const int numTasks =
          tbb::task_scheduler_init::default_num_threads();
        const size_t cacheSizeHint = 0;
        const size_t sizeOfScalar = sizeof(Scalar);

        RCP<ParameterList> params = parameterList ("NodeTsqr");
        params->set ("Num Tasks", numTasks,
                     "Number of tasks to use in the intranode parallel part "
                     "TSQR.  There is little/no performance penalty for mild "
                     "oversubscription, but a potential performance penalty "
                     "for undersubscription.");
        params->set ("Cache Size Hint", cacheSizeHint,
                    "Cache size hint in bytes (as a size_t) to use for "
                    "intranode TSQR.  If zero, TSQR will pick a reasonable "
                    "default.  See the documentation of SequentialTsqr for "
                     "a discussion of how to tune this parameter.");
        params->set ("Size of Scalar", sizeOfScalar);

        return params;
      }

      void
      setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
      {
        seq_.setParameterList (plist);

        if (! plist.is_null()) {
          const int defaultNumCores = 1; // A reasonable safe default value.
          int numTasks = plist->get ("Num Tasks", defaultNumCores);
          if (numTasks < 1) { // Default is no parallelism.
            plist->set ("Num Tasks", defaultNumCores);
          }
          numTasks_ = numTasks;
        }
      }

      /// \brief Number of tasks that TSQR will use to solve the problem.
      ///
      /// This is the number of subproblems into which to divide the
      /// main problem, in order to solve it in parallel.
      size_t ntasks() const { return numTasks_; }

      /// \brief Cache size hint (in bytes) used for the factorization.
      ///
      /// This may be different from the corresponding constructor
      /// argument, because TSQR may revise unreasonable suggestions
      /// into reasonable values.
      size_t cache_size_hint() const { return seq_.cache_size_hint(); }

      //! Fastest time over all tasks of the last SequentialTsqr::factor() call.
      double
      min_seq_factor_timing () const { return min_seq_factor_timing_; }
      //! Slowest time over all tasks of the last SequentialTsqr::factor() call.
      double
      max_seq_factor_timing () const { return max_seq_factor_timing_; }
      //! Fastest time over all tasks of the last SequentialTsqr::apply() call.
      double
      min_seq_apply_timing () const { return min_seq_apply_timing_; }
      //! Slowest time over all tasks of the last SequentialTsqr::apply() call.
      double
      max_seq_apply_timing () const { return max_seq_apply_timing_; }

      FactorOutput
      factor (const LocalOrdinal nrows,
              const LocalOrdinal ncols,
              Scalar A[],
              const LocalOrdinal lda,
              Scalar R[],
              const LocalOrdinal ldr,
              const bool contiguous_cache_blocks) const
      {
        using tbb::task;

        mat_view_type A_view (nrows, ncols, A, lda);
        // A_top will be modified in place by exactly one task, to
        // indicate the partition from which we may extract the R
        // factor after finishing the factorization.
        mat_view_type A_top;

        std::vector<SeqOutput> seq_output (ntasks());
        ParOutput par_output (ntasks(), std::vector<Scalar>(ncols));
        if (ntasks() < 1)
          {
            if (! A_view.empty())
              throw std::logic_error("Zero subproblems, but A not empty!");
            else // Return empty results
              return std::make_pair (seq_output, par_output);
          }

        double my_seq_timing = double(0);
        double min_seq_timing = double(0);
        double max_seq_timing = double(0);
        try {
          typedef FactorTask<LocalOrdinal, Scalar, TimerType> factor_task_t;

          // When the root task completes, A_top will be set to the
          // topmost partition of A.  We can then extract the R factor
          // from A_top.
          factor_task_t& root_task = *new( task::allocate_root() )
            factor_task_t(0, ntasks()-1, A_view, &A_top, seq_output,
                          par_output, seq_, my_seq_timing, min_seq_timing,
                          max_seq_timing, contiguous_cache_blocks);
          task::spawn_root_and_wait (root_task);
        } catch (tbb::captured_exception& ex) {
          // TBB can't guarantee on all systems that an exception
          // thrown in another thread will have its type correctly
          // propagated to this thread.  If it can't, then it captures
          // the exception as a tbb:captured_exception, and propagates
          // it to here.  It may be able to propagate the exception,
          // though, so be prepared for that.  We deal with the latter
          // case by allowing the exception to propagate.
          std::ostringstream os;
          os << "Intel TBB caught an exception, while computing the QR factor"
            "ization of a matrix A.  Unfortunately, its type information was "
            "lost, because the exception was thrown in another thread.  Its "
            "\"what()\" function returns the following string: " << ex.what();
          throw std::runtime_error (os.str());
        }

        // Copy the R factor out of A_top into R.
        seq_.extract_R (A_top.nrows(), A_top.ncols(), A_top.get(),
                        A_top.lda(), R, ldr, contiguous_cache_blocks);

        // Save the timings for future reference
        if (min_seq_timing < min_seq_factor_timing_)
          min_seq_factor_timing_ = min_seq_timing;
        if (max_seq_timing > max_seq_factor_timing_)
          max_seq_factor_timing_ = max_seq_timing;

        return std::make_pair (seq_output, par_output);
      }

      void
      apply (const ApplyType& apply_type,
             const LocalOrdinal nrows,
             const LocalOrdinal ncols_Q,
             const Scalar Q[],
             const LocalOrdinal ldq,
             const FactorOutput& factor_output,
             const LocalOrdinal ncols_C,
             Scalar C[],
             const LocalOrdinal ldc,
             const bool contiguous_cache_blocks) const
      {
        using tbb::task;

        if (apply_type.transposed())
          throw std::logic_error ("Applying Q^T and Q^H not implemented");

        const_mat_view_type Q_view (nrows, ncols_Q, Q, ldq);
        mat_view_type C_view (nrows, ncols_C, C, ldc);
        if (! apply_type.transposed())
          {
            array_top_blocks_t top_blocks (ntasks());
            build_partition_array (0, ntasks()-1, top_blocks, Q_view,
                                   C_view, contiguous_cache_blocks);
            double my_seq_timing = 0.0;
            double min_seq_timing = 0.0;
            double max_seq_timing = 0.0;
            try {
              typedef ApplyTask<LocalOrdinal, Scalar, TimerType> apply_task_t;
              apply_task_t& root_task =
                *new( task::allocate_root() )
                apply_task_t (0, ntasks()-1, Q_view, C_view, top_blocks,
                              factor_output, seq_, my_seq_timing,
                              min_seq_timing, max_seq_timing,
                              contiguous_cache_blocks);
              task::spawn_root_and_wait (root_task);
            } catch (tbb::captured_exception& ex) {
              std::ostringstream os;
              os << "Intel TBB caught an exception, while applying a Q factor "
                "computed previously by factor() to the matrix C.  Unfortunate"
                "ly, its type information was lost, because the exception was "
                "thrown in another thread.  Its \"what()\" function returns th"
                "e following string: " << ex.what();
              throw std::runtime_error (os.str());
            }

            // Save the timings for future reference
            if (min_seq_timing < min_seq_apply_timing_)
              min_seq_apply_timing_ = min_seq_timing;
            if (max_seq_timing > max_seq_apply_timing_)
              max_seq_apply_timing_ = max_seq_timing;
          }
      }


      void
      explicit_Q (const LocalOrdinal nrows,
                  const LocalOrdinal ncols_Q_in,
                  const Scalar Q_in[],
                  const LocalOrdinal ldq_in,
                  const FactorOutput& factor_output,
                  const LocalOrdinal ncols_Q_out,
                  Scalar Q_out[],
                  const LocalOrdinal ldq_out,
                  const bool contiguous_cache_blocks) const
      {
        using tbb::task;

        mat_view_type Q_out_view (nrows, ncols_Q_out, Q_out, ldq_out);
        try {
          typedef ExplicitQTask< LocalOrdinal, Scalar > explicit_Q_task_t;
          explicit_Q_task_t& root_task = *new( task::allocate_root() )
            explicit_Q_task_t (0, ntasks()-1, Q_out_view, seq_,
                               contiguous_cache_blocks);
          task::spawn_root_and_wait (root_task);
        } catch (tbb::captured_exception& ex) {
          std::ostringstream os;
          os << "Intel TBB caught an exception, while preparing to compute"
            " the explicit Q factor from a QR factorization computed previ"
            "ously by factor().  Unfortunately, its type information was l"
            "ost, because the exception was thrown in another thread.  Its"
            " \"what()\" function returns the following string: "
             << ex.what();
          throw std::runtime_error (os.str());
        }
        apply (ApplyType::NoTranspose,
               nrows, ncols_Q_in, Q_in, ldq_in, factor_output,
               ncols_Q_out, Q_out, ldq_out,
               contiguous_cache_blocks);
      }

      /// \brief Compute Q*B
      ///
      /// Compute matrix-matrix product Q*B, where Q is nrows by ncols
      /// and B is ncols by ncols.  Respect cache blocks of Q.
      void
      Q_times_B (const LocalOrdinal nrows,
                 const LocalOrdinal ncols,
                 Scalar Q[],
                 const LocalOrdinal ldq,
                 const Scalar B[],
                 const LocalOrdinal ldb,
                 const bool contiguous_cache_blocks) const
      {
        // Compute Q := Q*B in parallel.  This works much like
        // cache_block() (which see), in that each thread's instance
        // does not need to communicate with the others.
        try {
          using tbb::task;
          typedef RevealRankTask<LocalOrdinal, Scalar> rrtask_type;

          mat_view_type Q_view (nrows, ncols, Q, ldq);
          const_mat_view_type B_view (ncols, ncols, B, ldb);

          rrtask_type& root_task = *new( task::allocate_root() )
            rrtask_type (0, ntasks()-1, Q_view, B_view, seq_,
                         contiguous_cache_blocks);
          task::spawn_root_and_wait (root_task);
        } catch (tbb::captured_exception& ex) {
          std::ostringstream os;
          os << "Intel TBB caught an exception, while computing Q := Q*U.  "
            "Unfortunately, its type information was lost, because the "
            "exception was thrown in another thread.  Its \"what()\" function "
            "returns the following string: " << ex.what();
          throw std::runtime_error (os.str());
        }
      }


      /// Compute SVD \f$R = U \Sigma V^*\f$, not in place.  Use the
      /// resulting singular values to compute the numerical rank of R,
      /// with respect to the relative tolerance tol.  If R is full
      /// rank, return without modifying R.  If R is not full rank,
      /// overwrite R with \f$\Sigma \cdot V^*\f$.
      ///
      /// \return Numerical rank of R: 0 <= rank <= ncols.
      LocalOrdinal
      reveal_R_rank (const LocalOrdinal ncols,
                     Scalar R[],
                     const LocalOrdinal ldr,
                     Scalar U[],
                     const LocalOrdinal ldu,
                     const magnitude_type tol) const
      {
        return seq_.reveal_R_rank (ncols, R, ldr, U, ldu, tol);
      }

      /// \brief Rank-revealing decomposition
      ///
      /// Using the R factor from factor() and the explicit Q factor
      /// from explicit_Q(), compute the SVD of R (\f$R = U \Sigma
      /// V^*\f$).  R.  If R is full rank (with respect to the given
      /// relative tolerance tol), don't change Q or R.  Otherwise,
      /// compute \f$Q := Q \cdot U\f$ and \f$R := \Sigma V^*\f$ in
      /// place (the latter may be no longer upper triangular).
      ///
      /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq ncols\f$.
      ///
      LocalOrdinal
      reveal_rank (const LocalOrdinal nrows,
                   const LocalOrdinal ncols,
                   Scalar Q[],
                   const LocalOrdinal ldq,
                   Scalar R[],
                   const LocalOrdinal ldr,
                   const magnitude_type tol,
                   const bool contiguous_cache_blocks = false) const
      {
        // Take the easy exit if available.
        if (ncols == 0)
          return 0;

        Matrix<LocalOrdinal, Scalar> U (ncols, ncols, Scalar(0));
        const LocalOrdinal rank =
          reveal_R_rank (ncols, R, ldr, U.get(), U.ldu(), tol);

        if (rank < ncols)
          {
            // If R is not full rank: reveal_R_rank() already computed
            // the SVD \f$R = U \Sigma V^*\f$ of (the input) R, and
            // overwrote R with \f$\Sigma V^*\f$.  Now, we compute \f$Q
            // := Q \cdot U\f$, respecting cache blocks of Q.
            Q_times_B (nrows, ncols, Q, ldq, U.get(), U.lda(),
                       contiguous_cache_blocks);
          }
        return rank;
      }

      void
      cache_block (const LocalOrdinal nrows,
                   const LocalOrdinal ncols,
                   Scalar A_out[],
                   const Scalar A_in[],
                   const LocalOrdinal lda_in) const
      {
        using tbb::task;

        const_mat_view_type A_in_view (nrows, ncols, A_in, lda_in);
        // A_out won't have leading dimension lda_in, but that's OK,
        // as long as all the routines are told that A_out is
        // cache-blocked.
        mat_view_type A_out_view (nrows, ncols, A_out, lda_in);
        try {
          typedef CacheBlockTask< LocalOrdinal, Scalar > cache_block_task_t;
          cache_block_task_t& root_task = *new( task::allocate_root() )
            cache_block_task_t (0, ntasks()-1, A_out_view, A_in_view, seq_);
          task::spawn_root_and_wait (root_task);
        } catch (tbb::captured_exception& ex) {
          std::ostringstream os;
          os << "Intel TBB caught an exception, while cache-blocking a mat"
            "rix.  Unfortunately, its type information was lost, because t"
            "he exception was thrown in another thread.  Its \"what()\" fu"
            "nction returns the following string: " << ex.what();
          throw std::runtime_error (os.str());
        }
      }

      void
      un_cache_block (const LocalOrdinal nrows,
                      const LocalOrdinal ncols,
                      Scalar A_out[],
                      const LocalOrdinal lda_out,
                      const Scalar A_in[]) const
      {
        using tbb::task;

        // A_in doesn't have leading dimension lda_out, but that's OK,
        // as long as all the routines are told that A_in is cache-
        // blocked.
        const_mat_view_type A_in_view (nrows, ncols, A_in, lda_out);
        mat_view_type A_out_view (nrows, ncols, A_out, lda_out);
        try {
          typedef UnCacheBlockTask< LocalOrdinal, Scalar > un_cache_block_task_t;
          un_cache_block_task_t& root_task = *new( task::allocate_root() )
            un_cache_block_task_t (0, ntasks()-1, A_out_view, A_in_view, seq_);
          task::spawn_root_and_wait (root_task);
        } catch (tbb::captured_exception& ex) {
          std::ostringstream os;
          os << "Intel TBB caught an exception, while un-cache-blocking a "
            "matrix.  Unfortunately, its type information was lost, becaus"
            "e the exception was thrown in another thread.  Its \"what()\""
            " function returns the following string: " << ex.what();
          throw std::runtime_error (os.str());
        }
      }

      template< class MatrixViewType >
      MatrixViewType
      top_block (const MatrixViewType& C,
                 const bool contiguous_cache_blocks = false) const
      {
        return top_block_helper (0, ntasks()-1, C, contiguous_cache_blocks);
      }

      void
      fill_with_zeros (const LocalOrdinal nrows,
                       const LocalOrdinal ncols,
                       Scalar C[],
                       const LocalOrdinal ldc,
                       const bool contiguous_cache_blocks) const
      {
        using tbb::task;
        mat_view_type C_view (nrows, ncols, C, ldc);

        try {
          typedef FillWithZerosTask< LocalOrdinal, Scalar > fill_task_t;
          fill_task_t& root_task = *new( task::allocate_root() )
            fill_task_t (0, ntasks()-1, C_view, seq_, contiguous_cache_blocks);
          task::spawn_root_and_wait (root_task);
        } catch (tbb::captured_exception& ex) {
          std::ostringstream os;
          os << "Intel TBB caught an exception, while un-cache-blocking a "
            "matrix.  Unfortunately, its type information was lost, becaus"
            "e the exception was thrown in another thread.  Its \"what()\""
            " function returns the following string: " << ex.what();
          throw std::runtime_error (os.str());
        }
      }

    private:
      size_t numTasks_;
      TSQR::SequentialTsqr<LocalOrdinal, Scalar> seq_;
      TSQR::Combine<LocalOrdinal, Scalar> combine_;
      Partitioner<LocalOrdinal, Scalar> partitioner_;

      mutable double min_seq_factor_timing_;
      mutable double max_seq_factor_timing_;
      mutable double min_seq_apply_timing_;
      mutable double max_seq_apply_timing_;

      void
      build_partition_array (const size_t P_first,
                             const size_t P_last,
                             array_top_blocks_t& top_blocks,
                             const_mat_view_type& Q,
                             mat_view_type& C,
                             const bool contiguous_cache_blocks = false) const
      {
        if (P_first > P_last) {
          return;
        }
        else if (P_first == P_last) {
          const_mat_view_type Q_top = seq_.top_block (Q, contiguous_cache_blocks);
          mat_view_type C_top = seq_.top_block (C, contiguous_cache_blocks);
          top_blocks[P_first] =
            std::make_pair (const_mat_view_type (Q_top.ncols(), Q_top.ncols(),
                                                 Q_top.get(), Q_top.lda()),
                            mat_view_type (C_top.ncols(), C_top.ncols(),
                                           C_top.get(), C_top.lda()));
        }
        else {
          // Recurse on two intervals: [P_first, P_mid] and [P_mid+1, P_last]
          const size_t P_mid = (P_first + P_last) / 2;
          const_split_t Q_split =
            partitioner_.split (Q, P_first, P_mid, P_last,
                                contiguous_cache_blocks);
          split_t C_split =
            partitioner_.split (C, P_first, P_mid, P_last,
                                contiguous_cache_blocks);
          // The partitioner may decide that the current blocks Q
          // and C have too few rows to be worth splitting.  (The
          // partitioner should split both Q and C in the same way.)
          // In that case, Q_split.first should be the same block as
          // Q, and Q_split.second (the bottom block) will be empty.
          // Ditto for C_split.  We deal with this in the same way
          // as the base case (P_first == P_last) above.
          if (Q_split.second.empty() || Q_split.second.nrows() == 0) {
            const_mat_view_type Q_top =
              seq_.top_block (Q, contiguous_cache_blocks);
            mat_view_type C_top = seq_.top_block (C, contiguous_cache_blocks);
            top_blocks[P_first] =
              std::make_pair (const_mat_view_type (Q_top.ncols(), Q_top.ncols(),
                                                   Q_top.get(), Q_top.lda()),
                              mat_view_type (C_top.ncols(), C_top.ncols(),
                                             C_top.get(), C_top.lda()));
          }
          else {
            build_partition_array (P_first, P_mid, top_blocks,
                                   Q_split.first, C_split.first,
                                   contiguous_cache_blocks);
            build_partition_array (P_mid+1, P_last, top_blocks,
                                   Q_split.second, C_split.second,
                                   contiguous_cache_blocks);
          }
        }
      }
    };
  } // namespace TBB
} // namespace TSQR

#endif // __TSQR_TBB_TbbParallelTsqr_hpp
