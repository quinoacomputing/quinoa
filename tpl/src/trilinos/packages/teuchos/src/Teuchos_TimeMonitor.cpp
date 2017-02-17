// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_TableColumn.hpp"
#include "Teuchos_TableFormat.hpp"
#include <functional>
#ifdef HAVE_TEUCHOS_YAML_CPP
#  include <yaml-cpp/yaml.h>
#endif // HAVE_TEUCHOS_YAML_CPP


namespace Teuchos {
  /**
   * \class MaxLoc
   * \brief Teuchos version of MPI_MAXLOC.
   * \author Mark Hoemmen
   *
   * \tparam Ordinal The template parameter of \c Comm.
   *
   * \tparam Packet A type with value semantics; the type on which to
   *   reduce.
   *
   * \c MPI_MAXLOC is a standard reduction operator provided by the
   * MPI standard.  According to the standard, \c MPI_MAXLOC combines
   * the (value, index) pairs (u,i) and (v,j) into (w,j), where \f$w =
   * max(u,v)\f$, and
   * \f[
   *   k = \begin{cases}
   *     i         & \text{if $u > v$}, \\
   *     \min(i,j) & \text{if $u = v$}, \\
   *     j         & \text{if $u < v$}. \\
   * \end{cases}
   * \f]
   * This class implements the \c MPI_MAXLOC reduction operator for
   * the Teuchos communication wrappers.
   *
   * What happens to NaN ("Not a Number")?  A NaN is neither less
   * than, greater than, or equal to any floating-point number or any
   * NaN.  We can alter the above definition slightly so that a \c
   * MaxLoc reduction has a well-defined result in case the array
   * contains a NaN:
   *
   * \f[
   *   w = \begin{cases}
   *     u     & \text{if $u > v$}, \\
   *     v     & \text{if $u < v$}. \\
   *     u     & \text{otherwise}. \\
   *   \end{cases}
   * \f]
   * and
   * \f[
   *   k = \begin{cases}
   *     i         & \text{if $u > v$}, \\
   *     j         & \text{if $u < v$}. \\
   *     \min(i,j) & \text{otherwise}. \\
   *   \end{cases}
   * \f]
   * Defining \c MaxLoc in this way ensures that for any array
   * containing a NaN, the value (w) returned is the first NaN, and
   * the index (k) returned is the index of the first NaN.
   */
  template<class Ordinal, class ScalarType, class IndexType>
  class MaxLoc :
    public ValueTypeReductionOp<Ordinal, std::pair<ScalarType, IndexType> > {
  public:
    void
    reduce (const Ordinal count,
            const std::pair<ScalarType, IndexType> inBuffer[],
            std::pair<ScalarType, IndexType> inoutBuffer[]) const;
  };

  template<class Ordinal>
  class MaxLoc<Ordinal, double, int> :
    public ValueTypeReductionOp<Ordinal, std::pair<double, int> > {
  public:
    void
    reduce (const Ordinal count,
            const std::pair<double, int> inBuffer[],
            std::pair<double, int> inoutBuffer[]) const
    {
      for (Ordinal ind = 0; ind < count; ++ind) {
        const std::pair<double, int>& in = inBuffer[ind];
        std::pair<double, int>& inout = inoutBuffer[ind];

        if (in.first > inout.first) {
          inout.first = in.first;
          inout.second = in.second;
        } else if (in.first < inout.first) {
          // Don't need to do anything; inout has the values.
        } else { // equal, or at least one is NaN.
          inout.first = in.first;
          inout.second = std::min (in.second, inout.second);
        }
      }
    }
  };

  /** \class MinLoc
   * \brief Teuchos version of MPI_MINLOC.
   * \author Mark Hoemmen
   *
   * \tparam Ordinal The template parameter of \c Comm.
   *
   * \tparam Packet A type with value semantics; the type on which to
   *   reduce.
   *
   * \c MPI_MINLOC is a standard reduction operator provided by the
   * MPI standard.  According to the standard, \c MPI_MINLOC combines
   * the (value, index) pairs (u,i) and (v,j) into (w,j), where \f$w =
   * min(u,v)\f$, and
   * \f[
   *   k = \begin{cases}
   *     i         & \text{if $u < v$}, \\
   *     \min(i,j) & \text{if $u = v$}, \\
   *     j         & \text{if $u > v$}. \\
   *   \end{cases}
   * \f]
   * This class implements the \c MPI_MINLOC reduction operator for
   * the Teuchos communication wrappers.
   *
   * Refer to the note in the documentation of \c MaxLoc that
   * explains how we adjust the above definition to produce
   * well-defined results even if the array contains a NaN.
   */
  template<class Ordinal, class ScalarType, class IndexType>
  class MinLoc :
    public ValueTypeReductionOp<Ordinal, std::pair<ScalarType, IndexType> > {
  public:
    void
    reduce (const Ordinal count,
            const std::pair<ScalarType, IndexType> inBuffer[],
            std::pair<ScalarType, IndexType> inoutBuffer[]) const;
  };

  template<class Ordinal>
  class MinLoc<Ordinal, double, int> :
    public ValueTypeReductionOp<Ordinal, std::pair<double, int> > {
  public:
    void
    reduce (const Ordinal count,
            const std::pair<double, int> inBuffer[],
            std::pair<double, int> inoutBuffer[]) const
    {
      for (Ordinal ind = 0; ind < count; ++ind) {
        const std::pair<double, int>& in = inBuffer[ind];
        std::pair<double, int>& inout = inoutBuffer[ind];

        if (in.first < inout.first) {
          inout.first = in.first;
          inout.second = in.second;
        } else if (in.first > inout.first) {
          // Don't need to do anything; inout has the values.
        } else { // equal, or at least one is NaN.
          inout.first = in.first;
          inout.second = std::min (in.second, inout.second);
        }
      }
    }
  };

  // Typedef used internally by TimeMonitor::summarize() and its
  // helper functions.  The map is keyed on timer label (a string).
  // Each value is a pair: (total number of seconds over all calls to
  // that timer, total number of calls to that timer).
  typedef std::map<std::string, std::pair<double, int> > timer_map_t;

  TimeMonitor::TimeMonitor (Time& timer, bool reset)
    : PerformanceMonitorBase<Time>(timer, reset)
  {
    if (!isRecursiveCall()) counter().start(reset);
  }

  TimeMonitor::~TimeMonitor() {
    if (!isRecursiveCall()) counter().stop();
  }

  void
  TimeMonitor::zeroOutTimers()
  {
    const Array<RCP<Time> > timers = counters();

    // In debug mode, loop first to check whether any of the timers
    // are running, before resetting them.  This ensures that this
    // method satisfies the strong exception guarantee (either it
    // completes normally, or there are no side effects).
#ifdef TEUCHOS_DEBUG
    typedef Array<RCP<Time> >::size_type size_type;
    const size_type numTimers = timers.size();
    for (size_type i = 0; i < numTimers; ++i) {
      Time &timer = *timers[i];
      // We throw a runtime_error rather than a logic_error, because
      // logic_error suggests a bug in the implementation of
      // TimeMonitor.  Calling zeroOutTimers() when a timer is
      // running is not TimeMonitor's fault.
      TEUCHOS_TEST_FOR_EXCEPTION(timer.isRunning(), std::runtime_error,
                                 "The timer i = " << i << " with name \""
                                 << timer.name() << "\" is currently running and may not "
                                 "be reset.");
    }
#endif // TEUCHOS_DEBUG

    for (Array<RCP<Time> >::const_iterator it = timers.begin();
         it != timers.end(); ++it) {
      (*it)->reset ();
    }
  }

  // An anonymous namespace is the standard way of limiting linkage of
  // its contained routines to file scope.
  namespace {
    // \brief Return an "empty" local timer datum.
    //
    // "Empty" means the datum has zero elapsed time and zero call
    // count.  This function does not actually create a timer.
    //
    // \param name The timer's name.
    std::pair<std::string, std::pair<double, int> >
    makeEmptyTimerDatum (const std::string& name)
    {
      return std::make_pair (name, std::make_pair (double(0), int(0)));
    }

    // \fn collectLocalTimerData
    // \brief Collect and sort local timer data by timer names.
    //
    // \param localData [out] Map whose keys are the timer names, and
    //   whose value for each key is the total elapsed time (in
    //   seconds) and the call count for the timer with that name.
    //
    // \param localCounters [in] Timers from which to extract data.
    //
    // Extract the total elapsed time and call count from each timer
    // in the given array.  Merge results for timers with duplicate
    // labels, by summing their total elapsed times and call counts
    // pairwise.
    void
    collectLocalTimerData (timer_map_t& localData,
                           ArrayView<const RCP<Time> > localCounters)
    {
      using std::make_pair;
      typedef timer_map_t::const_iterator const_iter_t;
      typedef timer_map_t::iterator iter_t;

      timer_map_t theLocalData;
      for (ArrayView<const RCP<Time> >::const_iterator it = localCounters.begin();
           it != localCounters.end(); ++it) {
        const std::string& name = (*it)->name();
        const double timing = (*it)->totalElapsedTime();
        const int numCalls = (*it)->numCalls();

        // Merge timers with duplicate labels, by summing their
        // total elapsed times and call counts.
        iter_t loc = theLocalData.find (name);
        if (loc == theLocalData.end()) {
          // Use loc as an insertion location hint.
          theLocalData.insert (loc, make_pair (name, make_pair (timing, numCalls)));
        }
        else {
          loc->second.first += timing;
          loc->second.second += numCalls;
        }
      }
      // This avoids copying the map, and also makes this method
      // satisfy the strong exception guarantee.
      localData.swap (theLocalData);
    }

    // \brief Locally filter out timer data with zero call counts.
    //
    // \param timerData [in/out]
    void
    filterZeroData (timer_map_t& timerData)
    {
      timer_map_t newTimerData;
      for (timer_map_t::const_iterator it = timerData.begin();
           it != timerData.end(); ++it) {
        if (it->second.second > 0) {
          newTimerData[it->first] = it->second;
        }
      }
      timerData.swap (newTimerData);
    }

    /// \fn collectLocalTimerDataAndNames
    /// \brief Collect the local timer data and timer names from the timers.
    ///
    /// \param localTimerData [out] Timer data extracted from
    ///   localTimers.  Contents on input are ignored and overwritten.
    ///   See the documentation of \c collectLocalTimerData().
    ///
    /// \param localTimerNames [out] On output: names of timers
    ///   extracted from \c localTimers, in the same order as the keys
    ///   of \c localTimerData.  Resized as necessary.  Contents on
    ///   input are ignored and overwritten.
    ///
    /// \param localTimers [in] (Local) timers from which to extract data.
    ///
    /// \param writeZeroTimers [in] If true, do not include timers
    ///   with zero call counts in the \c localTimerData and \c
    ///   localTimerNames output.
    void
    collectLocalTimerDataAndNames (timer_map_t& localTimerData,
                                   Array<std::string>& localTimerNames,
                                   ArrayView<const RCP<Time> > localTimers,
                                   const bool writeZeroTimers)
    {
      // Collect and sort local timer data by timer names.
      collectLocalTimerData (localTimerData, localTimers);

      // Filter out zero data locally first.  This ensures that if we
      // are writing global stats, and if a timer name exists in the
      // set of global names, then that timer has a nonzero call count
      // on at least one MPI process.
      if (! writeZeroTimers) {
        filterZeroData (localTimerData);
      }

      // Extract the set of local timer names.  The std::map keeps
      // them sorted alphabetically.
      localTimerNames.reserve (localTimerData.size());
      for (timer_map_t::const_iterator it = localTimerData.begin();
           it != localTimerData.end(); ++it) {
        localTimerNames.push_back (it->first);
      }
    }

    /// \brief Merge local timer data into global data.
    ///
    /// Call this method in \c summarize() only if the \c
    /// writeGlobalStats argument is true.
    ///
    /// \param globalTimerData [out] Result of merging localTimerData
    ///   over the processes in the given communicator.
    ///
    /// \param globalTimerNames [out] Names of the timers in
    ///   globalTimerData; same as the keys of the map.
    ///
    /// \param localTimerData [in/out] On input: the first return
    ///   value of \c collectLocalTimerDataAndNames().  On output, if
    ///   writeZeroTimers is true, data for timers with zero call
    ///   counts will be removed.
    ///
    /// \param localTimerNames [in/out] On input: the second return
    ///   value of \c collectLocalTimerDataAndNames().  On output, if
    ///   writeZeroTimers is true, names of timers with zero call
    ///   counts will be removed.
    ///
    /// \param comm [in] Communicator over which to merge.
    ///
    /// \param alwaysWriteLocal [in] If true, and if the local set of
    ///   timers differs from the global set of timers (either the
    ///   union or the intersection, depending on \c setOp), Proc 0
    ///   will create corresponding timer data in localTimerData with
    ///   zero elapsed times and call counts.  (This pads the output
    ///   so it fits in a tabular form.)
    ///
    /// \param setOp [in] If \c Intersection, compute the intersection
    ///   of all created timers over all processes in the
    ///   communicator.  If \c Union, compute the union of all created
    ///   timers over all processes in the communicator.
    void
    collectGlobalTimerData (timer_map_t& globalTimerData,
                            Array<std::string>& globalTimerNames,
                            timer_map_t& localTimerData,
                            Array<std::string>& localTimerNames,
                            Ptr<const Comm<int> > comm,
                            const bool alwaysWriteLocal,
                            const ECounterSetOp setOp)
    {
      // There may be some global timers that are not local timers on
      // the calling MPI process(es).  In that case, if
      // alwaysWriteLocal is true, then we need to fill in the
      // "missing" local timers.  That will ensure that both global
      // and local timer columns in the output table have the same
      // number of rows.  The collectLocalTimerDataAndNames() method
      // may have already filtered out local timers with zero call
      // counts (if its writeZeroTimers argument was false), but we
      // won't be filtering again.  Thus, any local timer data we
      // insert here won't get filtered out.
      //
      // Note that calling summarize() with writeZeroTimers == false
      // will still do what it says, even if we insert local timers
      // with zero call counts here.

      // This does the correct and inexpensive thing (just copies the
      // timer data) if numProcs == 1.  Otherwise, it initiates a
      // communication with \f$O(\log P)\f$ messages along the
      // critical path, where \f$P\f$ is the number of participating
      // processes.
      mergeCounterNames (*comm, localTimerNames, globalTimerNames, setOp);

#ifdef TEUCHOS_DEBUG
      {
        // Sanity check that all processes have the name number of
        // global timer names.
        const timer_map_t::size_type myNumGlobalNames = globalTimerNames.size();
        timer_map_t::size_type minNumGlobalNames = 0;
        timer_map_t::size_type maxNumGlobalNames = 0;
        reduceAll (*comm, REDUCE_MIN, myNumGlobalNames,
                   outArg (minNumGlobalNames));
        reduceAll (*comm, REDUCE_MAX, myNumGlobalNames,
                   outArg (maxNumGlobalNames));
        TEUCHOS_TEST_FOR_EXCEPTION(minNumGlobalNames != maxNumGlobalNames,
          std::logic_error, "Min # global timer names = " << minNumGlobalNames
          << " != max # global timer names = " << maxNumGlobalNames
          << ".  Please report this bug to the Teuchos developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(myNumGlobalNames != minNumGlobalNames,
          std::logic_error, "My # global timer names = " << myNumGlobalNames
          << " != min # global timer names = " << minNumGlobalNames
          << ".  Please report this bug to the Teuchos developers.");
      }
#endif // TEUCHOS_DEBUG

      // mergeCounterNames() just merges the counters' names, not
      // their actual data.  Now we need to fill globalTimerData with
      // this process' timer data for the timers in globalTimerNames.
      //
      // All processes need the full list of global timers, since
      // there may be some global timers that are not local timers.
      // That's why mergeCounterNames() has to be an all-reduce, not
      // just a reduction to Proc 0.
      //
      // Insertion optimization: if the iterator given to map::insert
      // points right before where we want to insert, insertion is
      // O(1).  globalTimerNames is sorted, so feeding the iterator
      // output of map::insert into the next invocation's input should
      // make the whole insertion O(N) where N is the number of
      // entries in globalTimerNames.
      timer_map_t::iterator globalMapIter = globalTimerData.begin();
      timer_map_t::iterator localMapIter;
      for (Array<string>::const_iterator it = globalTimerNames.begin();
           it != globalTimerNames.end(); ++it) {
        const std::string& globalName = *it;
        localMapIter = localTimerData.find (globalName);

        if (localMapIter == localTimerData.end()) {
          if (alwaysWriteLocal) {
            // If there are some global timers that are not local
            // timers, and if we want to print local timers, we insert
            // a local timer datum with zero elapsed time and zero
            // call count into localTimerData as well.  This will
            // ensure that both global and local timer columns in the
            // output table have the same number of rows.
            //
            // We really only need to do this on Proc 0, which is the
            // only process that currently may print local timers.
            // However, we do it on all processes, just in case
            // someone later wants to modify this function to print
            // out local timer data for some process other than Proc
            // 0.  This extra computation won't affect the cost along
            // the critical path, for future computations in which
            // Proc 0 participates.
            localMapIter = localTimerData.insert (localMapIter, makeEmptyTimerDatum (globalName));

            // Make sure the missing global name gets added to the
            // list of local names.  We'll re-sort the list of local
            // names below.
            localTimerNames.push_back (globalName);
          }
          // There's a global timer that's not a local timer.  Add it
          // to our pre-merge version of the global timer data so that
          // we can safely merge the global timer data later.
          globalMapIter = globalTimerData.insert (globalMapIter, makeEmptyTimerDatum (globalName));
        }
        else {
          // We have this global timer name in our local timer list.
          // Fill in our pre-merge version of the global timer data
          // with our local data.
          globalMapIter = globalTimerData.insert (globalMapIter, std::make_pair (globalName, localMapIter->second));
        }
      }

      if (alwaysWriteLocal) {
        // Re-sort the list of local timer names, since we may have
        // inserted "missing" names above.
        std::sort (localTimerNames.begin(), localTimerNames.end());
      }

#ifdef TEUCHOS_DEBUG
      {
        // Sanity check that all processes have the name number of
        // global timers.
        const timer_map_t::size_type myNumGlobalTimers = globalTimerData.size();
        timer_map_t::size_type minNumGlobalTimers = 0;
        timer_map_t::size_type maxNumGlobalTimers = 0;
        reduceAll (*comm, REDUCE_MIN, myNumGlobalTimers,
                   outArg (minNumGlobalTimers));
        reduceAll (*comm, REDUCE_MAX, myNumGlobalTimers,
                   outArg (maxNumGlobalTimers));
        TEUCHOS_TEST_FOR_EXCEPTION(minNumGlobalTimers != maxNumGlobalTimers,
                                   std::logic_error, "Min # global timers = " << minNumGlobalTimers
                                   << " != max # global timers = " << maxNumGlobalTimers
                                   << ".  Please report this bug to the Teuchos developers.");
        TEUCHOS_TEST_FOR_EXCEPTION(myNumGlobalTimers != minNumGlobalTimers,
                                   std::logic_error, "My # global timers = " << myNumGlobalTimers
                                   << " != min # global timers = " << minNumGlobalTimers
                                   << ".  Please report this bug to the Teuchos developers.");
      }
#endif // TEUCHOS_DEBUG
    }

    /// \brief Compute global timer statistics.
    ///
    /// Currently, this function computes the "MinOverProcs",
    /// "MeanOverProcs", "MaxOverProcs", and "MeanOverCallCounts"
    /// timings for each timer.  Along with the min / max timing comes
    /// the call count of the process who had the min / max.  (If more
    /// than one process had the min / max, then the call count on the
    /// process with the smallest rank is reported.)
    ///
    /// "MeanOverProcs" is the arithmetic mean of the cumulative
    /// timings over processes.  It ignores the call counts.  It
    /// includes the mean call count, which may be fractional, and has
    /// no particular connection to the mean timing.
    ///
    /// The "MeanOverCallCounts" is an arithmetic mean of timings that
    /// accounts for call counts.  Each timing is the sum over all
    /// calls.  Thus, this mean equals the sum of the timing over all
    /// processes, divided by the sum of the call counts over all
    /// processes for that timing.  (We compute it a bit differently
    /// to help prevent overflow.)  Along with the mean timing comes
    /// the same mean call count as mentioned above.
    ///
    /// \param statData [out] On output: Global timer statistics.  See
    ///   the \c stat_map_type typedef documentation for an explanation
    ///   of the data structure.
    ///
    /// \param statNames [out] On output: Each value in the statData
    ///   map is a vector.  That vector v has the same number of
    ///   entries as statNames.  statNames[k] is the name of the
    ///   statistic (e.g., "Min", "MeanOverProcs", "Max", or
    ///   "MeanOverCallCounts") stored as v[k].
    ///
    /// \param comm [in] Communicator over which to compute statistics.
    ///
    /// \param globalTimerData [in] Output with the same name of the
    ///   \c collectGlobalTimerData() function.  That function assures
    ///   that all processes have the same keys stored in this map.
    void
    computeGlobalTimerStats (stat_map_type& statData,
                             std::vector<std::string>& statNames,
                             Ptr<const Comm<int> > comm,
                             const timer_map_t& globalTimerData)
    {
      const int numTimers = static_cast<int> (globalTimerData.size());
      const int numProcs = comm->getSize();

      // Extract pre-reduction timings and call counts into a
      // sequential array.  This array will be in the same order as
      // the global timer names are in the map.
      Array<std::pair<double, int> > timingsAndCallCounts;
      timingsAndCallCounts.reserve (numTimers);
      for (timer_map_t::const_iterator it = globalTimerData.begin();
           it != globalTimerData.end(); ++it) {
        timingsAndCallCounts.push_back (it->second);
      }

      // For each timer name, compute the min timing and its
      // corresponding call count.  If two processes have the same
      // timing but different call counts, the minimum call count will
      // be used.
      Array<std::pair<double, int> > minTimingsAndCallCounts (numTimers);
      if (numTimers > 0) {
        reduceAll (*comm, MinLoc<int, double, int>(), numTimers,
                   &timingsAndCallCounts[0], &minTimingsAndCallCounts[0]);
      }

      // For each timer name, compute the max timing and its
      // corresponding call count.  If two processes have the same
      // timing but different call counts, the minimum call count will
      // be used.
      Array<std::pair<double, int> > maxTimingsAndCallCounts (numTimers);
      if (numTimers > 0) {
        reduceAll (*comm, MaxLoc<int, double, int>(), numTimers,
                   &timingsAndCallCounts[0], &maxTimingsAndCallCounts[0]);
      }

      // For each timer name, compute the mean-over-processes timing,
      // the mean call count, and the mean-over-call-counts timing.
      // The mean call count is reported as a double to allow a
      // fractional value.
      //
      // Each local timing is really the total timing over all local
      // invocations.  The number of local invocations is the call
      // count.  Thus, the mean-over-call-counts timing is the sum of
      // all the timings (over all processes), divided by the sum of
      // all the call counts (over all processes).  We compute it in a
      // different way to over unnecessary overflow.
      Array<double> meanOverCallCountsTimings (numTimers);
      Array<double> meanOverProcsTimings (numTimers);
      Array<double> meanCallCounts (numTimers);
      {
        // When summing, first scale by the number of processes.  This
        // avoids unnecessary overflow, and also gives us the mean
        // call count automatically.
        Array<double> scaledTimings (numTimers);
        Array<double> scaledCallCounts (numTimers);
        const double P = static_cast<double> (numProcs);
        for (int k = 0; k < numTimers; ++k) {
          const double timing = timingsAndCallCounts[k].first;
          const double callCount = static_cast<double> (timingsAndCallCounts[k].second);

          scaledTimings[k] = timing / P;
          scaledCallCounts[k] = callCount / P;
        }
        if (numTimers > 0) {
          reduceAll (*comm, REDUCE_SUM, numTimers, &scaledTimings[0],
                     &meanOverProcsTimings[0]);
          reduceAll (*comm, REDUCE_SUM, numTimers, &scaledCallCounts[0],
                     &meanCallCounts[0]);
        }
        // We don't have to undo the scaling for the mean timings;
        // just divide by the scaled call count.
        for (int k = 0; k < numTimers; ++k) {
          meanOverCallCountsTimings[k] = meanOverProcsTimings[k] / meanCallCounts[k];
        }
      }

      // Reformat the data into the map of statistics.  Be sure that
      // each value (the std::vector of (timing, call count) pairs,
      // each entry of which is a different statistic) preserves the
      // order of statNames.
      statNames.resize (4);
      statNames[0] = "MinOverProcs";
      statNames[1] = "MeanOverProcs";
      statNames[2] = "MaxOverProcs";
      statNames[3] = "MeanOverCallCounts";

      stat_map_type::iterator statIter = statData.end();
      timer_map_t::const_iterator it = globalTimerData.begin();
      for (int k = 0; it != globalTimerData.end(); ++k, ++it) {
        std::vector<std::pair<double, double> > curData (4);
        curData[0] = minTimingsAndCallCounts[k];
        curData[1] = std::make_pair (meanOverProcsTimings[k], meanCallCounts[k]);
        curData[2] = maxTimingsAndCallCounts[k];
        curData[3] = std::make_pair (meanOverCallCountsTimings[k], meanCallCounts[k]);

        // statIter gives an insertion location hint that makes each
        // insertion O(1), since we remember the location of the last
        // insertion.
        statIter = statData.insert (statIter, std::make_pair (it->first, curData));
      }
    }


    /// \brief Get a default communicator appropriate for the environment.
    ///
    /// If Trilinos was configured with MPI support, and if MPI has
    /// been initialized, return a wrapped MPI_COMM_WORLD.  If
    /// Trilinos was configured with MPI support, and if MPI has not
    /// yet been initialized, return a serial communicator (containing
    /// one process).  If Trilinos was <i>not</i> configured with MPI
    /// support, return a serial communicator.
    ///
    /// Rationale: Callers may or may not have initialized MPI before
    /// calling this method.  Just because they built with MPI,
    /// doesn't mean they want to use MPI.  It's not my responsibility
    /// to initialize MPI for them, and I don't have the context I
    /// need in order to do so anyway.  Thus, if Trilinos was built
    /// with MPI and MPI has not yet been initialized, this method
    /// returns a "serial" communicator.
    RCP<const Comm<int> >
    getDefaultComm ()
    {
      // The default communicator.  If Trilinos was built with MPI
      // enabled, this should be MPI_COMM_WORLD.  (If MPI has not yet
      // been initialized, it's not valid to use the communicator!)
      // Otherwise, this should be a "serial" (no MPI, one "process")
      // communicator.
      RCP<const Comm<int> > comm = DefaultComm<int>::getComm ();

#ifdef HAVE_MPI
      {
        int mpiHasBeenStarted = 0;
        MPI_Initialized (&mpiHasBeenStarted);
        if (! mpiHasBeenStarted) {
          // Make pComm a new "serial communicator."
          comm = rcp_implicit_cast<const Comm<int> > (rcp (new SerialComm<int> ()));
        }
      }
#endif // HAVE_MPI
      return comm;
    }

  } // namespace (anonymous)


  void
  TimeMonitor::computeGlobalTimerStatistics (stat_map_type& statData,
                                             std::vector<std::string>& statNames,
                                             Ptr<const Comm<int> > comm,
                                             const ECounterSetOp setOp)
  {
    // Collect local timer data and names.  Filter out timers with
    // zero call counts if writeZeroTimers is false.
    timer_map_t localTimerData;
    Array<std::string> localTimerNames;
    const bool writeZeroTimers = false;
    collectLocalTimerDataAndNames (localTimerData, localTimerNames,
                                   counters(), writeZeroTimers);
    // Merge the local timer data and names into global timer data and
    // names.
    timer_map_t globalTimerData;
    Array<std::string> globalTimerNames;
    const bool alwaysWriteLocal = false;
    collectGlobalTimerData (globalTimerData, globalTimerNames,
                            localTimerData, localTimerNames,
                            comm, alwaysWriteLocal, setOp);
    // Compute statistics on the data.
    computeGlobalTimerStats (statData, statNames, comm, globalTimerData);
  }


  void
  TimeMonitor::summarize (Ptr<const Comm<int> > comm,
                          std::ostream& out,
                          const bool alwaysWriteLocal,
                          const bool writeGlobalStats,
                          const bool writeZeroTimers,
                          const ECounterSetOp setOp)
  {
    //
    // We can't just call computeGlobalTimerStatistics(), since
    // summarize() has different options that affect whether global
    // statistics are computed and printed.
    //
    const int numProcs = comm->getSize();
    const int myRank = comm->getRank();

    // Collect local timer data and names.  Filter out timers with
    // zero call counts if writeZeroTimers is false.
    timer_map_t localTimerData;
    Array<std::string> localTimerNames;
    collectLocalTimerDataAndNames (localTimerData, localTimerNames,
                                   counters(), writeZeroTimers);

    // If we're computing global statistics, merge the local timer
    // data and names into global timer data and names, and compute
    // global timer statistics.  Otherwise, leave the global data
    // empty.
    timer_map_t globalTimerData;
    Array<std::string> globalTimerNames;
    stat_map_type statData;
    std::vector<std::string> statNames;
    if (writeGlobalStats) {
      collectGlobalTimerData (globalTimerData, globalTimerNames,
                              localTimerData, localTimerNames,
                              comm, alwaysWriteLocal, setOp);
      // Compute statistics on the data, but only if the communicator
      // contains more than one process.  Otherwise, statistics don't
      // make sense and we don't print them (see below).
      if (numProcs > 1) {
        computeGlobalTimerStats (statData, statNames, comm, globalTimerData);
      }
    }

    // Precision of floating-point numbers in the table.
    const int precision = format().precision();

    // All columns of the table, in order.
    Array<TableColumn> tableColumns;

    // Labels of all the columns of the table.
    // We will append to this when we add each column.
    Array<std::string> titles;

    // Widths (in number of characters) of each column.
    // We will append to this when we add each column.
    Array<int> columnWidths;

    // Table column containing all timer names.  If writeGlobalStats
    // is true, we use the global timer names, otherwise we use the
    // local timer names.  We build the table on all processes
    // redundantly, but only print on Rank 0.
    {
      titles.append ("Timer Name");

      // The column labels depend on whether we are computing global statistics.
      TableColumn nameCol (writeGlobalStats ? globalTimerNames : localTimerNames);
      tableColumns.append (nameCol);

      // Each column is as wide as it needs to be to hold both its
      // title and all of the column data.  This column's title is the
      // current last entry of the titles array.
      columnWidths.append (format().computeRequiredColumnWidth (titles.back(), nameCol));
    }

    // Table column containing local timer stats, if applicable.  We
    // only write local stats if asked, only on MPI Proc 0, and only
    // if there is more than one MPI process in the communicator
    // (otherwise local stats == global stats, so we just print the
    // global stats).  In this case, we've padded the local data on
    // Proc 0 if necessary to match the global timer list, so that the
    // columns have the same number of rows.
    if (alwaysWriteLocal && numProcs > 1 && myRank == 0) {
      titles.append ("Local time (num calls)");

      // Copy local timer data out of the array-of-structs into
      // separate arrays, for display in the table.
      Array<double> localTimings;
      Array<double> localNumCalls;
      for (timer_map_t::const_iterator it = localTimerData.begin();
           it != localTimerData.end(); ++it) {
        localTimings.push_back (it->second.first);
        localNumCalls.push_back (static_cast<double> (it->second.second));
      }
      TableColumn timeAndCalls (localTimings, localNumCalls, precision, true);
      tableColumns.append (timeAndCalls);
      columnWidths.append (format().computeRequiredColumnWidth (titles.back(), timeAndCalls));
    }

    if (writeGlobalStats) {
      // If there's only 1 process in the communicator, don't display
      // statistics; statistics don't make sense in that case.  Just
      // display the timings and call counts.  If there's more than 1
      // process, do display statistics.
      if (numProcs == 1) {
        // Extract timings and the call counts from globalTimerData.
        Array<double> globalTimings;
        Array<double> globalNumCalls;
        for (timer_map_t::const_iterator it = globalTimerData.begin();
             it != globalTimerData.end(); ++it) {
          globalTimings.push_back (it->second.first);
          globalNumCalls.push_back (static_cast<double> (it->second.second));
        }
        // Print the table column.
        titles.append ("Global time (num calls)");
        TableColumn timeAndCalls (globalTimings, globalNumCalls, precision, true);
        tableColumns.append (timeAndCalls);
        columnWidths.append (format().computeRequiredColumnWidth (titles.back(), timeAndCalls));
      }
      else { // numProcs > 1
        // Print a table column for each statistic.  statNames and
        // each value in statData use the same ordering, so we can
        // iterate over valid indices of statNames to display the
        // statistics in the right order.
        const timer_map_t::size_type numGlobalTimers = globalTimerData.size();
        for (std::vector<std::string>::size_type statInd = 0; statInd < statNames.size(); ++statInd) {
          // Extract lists of timings and their call counts for the
          // current statistic.
          Array<double> statTimings (numGlobalTimers);
          Array<double> statCallCounts (numGlobalTimers);
          stat_map_type::const_iterator it = statData.begin();
          for (int k = 0; it != statData.end(); ++it, ++k) {
            statTimings[k] = (it->second[statInd]).first;
            statCallCounts[k] = (it->second[statInd]).second;
          }
          // Print the table column.
          const std::string& statisticName = statNames[statInd];
          const std::string titleString = statisticName;
          titles.append (titleString);
          TableColumn timeAndCalls (statTimings, statCallCounts, precision, true);
          tableColumns.append (timeAndCalls);
          columnWidths.append (format().computeRequiredColumnWidth (titles.back(), timeAndCalls));
        }
      }
    }

    // Print the whole table to the given output stream on MPI Rank 0.
    format().setColumnWidths (columnWidths);
    if (myRank == 0) {
      std::ostringstream theTitle;
      theTitle << "TimeMonitor results over " << numProcs << " processor"
               << (numProcs > 1 ? "s" : "");
      format().writeWholeTable (out, theTitle.str(), titles, tableColumns);
    }
  }

  void
  TimeMonitor::summarize (std::ostream &out,
                          const bool alwaysWriteLocal,
                          const bool writeGlobalStats,
                          const bool writeZeroTimers,
                          const ECounterSetOp setOp)
  {
    // The default communicator.  If Trilinos was built with MPI
    // enabled, this should be MPI_COMM_WORLD.  Otherwise, this should
    // be a "serial" (no MPI, one "process") communicator.
    RCP<const Comm<int> > comm = getDefaultComm();

    summarize (comm.ptr(), out, alwaysWriteLocal,
               writeGlobalStats, writeZeroTimers, setOp);
  }

  void
  TimeMonitor::computeGlobalTimerStatistics (stat_map_type& statData,
                                             std::vector<std::string>& statNames,
                                             const ECounterSetOp setOp)
  {
    // The default communicator.  If Trilinos was built with MPI
    // enabled, this should be MPI_COMM_WORLD.  Otherwise, this should
    // be a "serial" (no MPI, one "process") communicator.
    RCP<const Comm<int> > comm = getDefaultComm();

    computeGlobalTimerStatistics (statData, statNames, comm.ptr(), setOp);
  }


  void TimeMonitor::
  summarizeToYaml (Ptr<const Comm<int> > comm, std::ostream &out)
  {
#ifdef HAVE_TEUCHOS_YAML_CPP
    // const bool writeGlobalStats = true;
    // const bool writeZeroTimers = true;
    // const bool alwaysWriteLocal = false;
    const ECounterSetOp setOp = Intersection;

    stat_map_type statData;
    std::vector<std::string> statNames;
    computeGlobalTimerStatistics (statData, statNames, setOp);

    const int numProcs = comm->getSize();
    const int myRank = comm->getRank();

    if (myRank == 0) {
      YAML::Emitter emi;
      emi << YAML::BeginDoc; // Begin YAML output
      emi << "Teuchos::TimeMonitor timing results";
      emi << YAML::BeginMap // Begin timing results map
          << YAML::Key << "Number of processes"
          << YAML::Value << numProcs
          << YAML::Key << "Global timer statistics"
          << YAML::Value;
      // For each timer name, print all its statistics.
      emi << YAML::BeginMap; // Begin timer names
      for (stat_map_type::const_iterator statDataIter = statData.begin();
           statDataIter != statData.end(); ++statDataIter) {
        // Key: Timer's name
        emi << YAML::Key << statDataIter->first;
        // Value: The timer's statistics, as a map.
        emi << YAML::Value << YAML::BeginMap; // Begin current timer's statistics
        for (std::vector<std::string>::size_type statInd = 0;
             statInd < statNames.size (); ++statInd) {
          // Key: current statistic's name.
          emi << YAML::Key << statNames[statInd]
              << YAML::Value;
          // Value is a map: "Time (s)" => time in seconds for current
          // statistic, "Call count" => call count for current statistic.
          const double curTime = (statDataIter->second)[statInd].first;
          const double curCallCount = (statDataIter->second)[statInd].second;
          emi << YAML::BeginMap
              << YAML::Key << "Time (s)"
              << YAML::Value << curTime
              << YAML::Key << "Call count"
              << YAML::Value << curCallCount
              << YAML::EndMap;
        }
        emi << YAML::EndMap; // End current timer's statistics
      }
      emi << YAML::EndMap; // End timer names
      emi << YAML::EndMap; // End timing results map
      emi << YAML::EndDoc; // End YAML output

      // Write YAML output to the given output stream.
      out << emi.c_str ();
    }
#else  // Don't HAVE_TEUCHOS_YAML_CPP
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Teuchos::TimeMonitor: YAML output currently requires building Trilinos with the yaml-cpp library.  Please download and install yaml-cpp from http://code.google.com/p/yaml-cpp/.  Then, enable yaml-cpp support when building Trilinos: 1. Set the CMake Boolean option TPL_ENABLE_yaml-cpp to ON.  2. Set the CMake option yaml-cpp_INCLUDE_DIRS to the path of the yaml-cpp include files (not including the yaml-cpp directory in include/).  3. Set the CMake option yaml-cpp_LIBRARY_DIRS to the location of the yaml-cpp library.  4. Clear the CMake cache if necesssary.  5. Run CMake again and rebuild Trilinos.");
#endif // HAVE_TEUCHOS_YAML_CPP
  }

  void TimeMonitor::
  summarizeToYaml (std::ostream &out)
  {
    // The default communicator.  If Trilinos was built with MPI
    // enabled, this should be MPI_COMM_WORLD.  Otherwise, this should
    // be a "serial" (no MPI, one "process") communicator.
    RCP<const Comm<int> > comm = getDefaultComm ();

    summarizeToYaml (comm.ptr (), out);
  }

  // Default value is false.  We'll set to true once
  // setReportParameters() completes successfully.
  bool TimeMonitor::setParams_ = false;

  // We have to declare all of these here in order to avoid linker errors.
  TimeMonitor::ETimeMonitorReportFormat TimeMonitor::reportFormat_ =
             TimeMonitor::REPORT_FORMAT_TABLE;
  ECounterSetOp TimeMonitor::setOp_ = Intersection;
  bool TimeMonitor::alwaysWriteLocal_ = false;
  bool TimeMonitor::writeGlobalStats_ = true;
  bool TimeMonitor::writeZeroTimers_ = true;

  void
  TimeMonitor::setReportFormatParameter (ParameterList& plist)
  {
    const std::string name ("Report format");
    const std::string defaultValue ("Table");
    const std::string docString ("Output format for report of timer statistics");
    Array<std::string> strings;
    Array<std::string> docs;
    Array<ETimeMonitorReportFormat> values;

    strings.push_back ("YAML");
    docs.push_back ("YAML (see yaml.org) format");
    values.push_back (REPORT_FORMAT_YAML);
    strings.push_back ("Table");
    docs.push_back ("Tabular format via Teuchos::TableFormat");
    values.push_back (REPORT_FORMAT_TABLE);

    setStringToIntegralParameter<ETimeMonitorReportFormat> (name, defaultValue,
                                                            docString,
                                                            strings (), docs (),
                                                            values (), &plist);
  }

  void
  TimeMonitor::setSetOpParameter (ParameterList& plist)
  {
    const std::string name ("How to merge timer sets");
    const std::string defaultValue ("Intersection");
    const std::string docString ("How to merge differing sets of timers "
                                 "across processes");
    Array<std::string> strings;
    Array<std::string> docs;
    Array<ECounterSetOp> values;

    strings.push_back ("Intersection");
    docs.push_back ("Compute intersection of timer sets over processes");
    values.push_back (Intersection);
    strings.push_back ("Union");
    docs.push_back ("Compute union of timer sets over processes");
    values.push_back (Union);

    setStringToIntegralParameter<ECounterSetOp> (name, defaultValue, docString,
                                                 strings (), docs (), values (),
                                                 &plist);
  }

  RCP<const ParameterList>
  TimeMonitor::getValidReportParameters ()
  {
    // Our implementation favors recomputation over persistent
    // storage.  That is, we simply recreate the list every time we
    // need it.
    RCP<ParameterList> plist = parameterList ("TimeMonitor::report");

    const bool alwaysWriteLocal = false;
    const bool writeGlobalStats = true;
    const bool writeZeroTimers = true;

    setReportFormatParameter (*plist);
    setSetOpParameter (*plist);
    plist->set ("alwaysWriteLocal", alwaysWriteLocal,
                "Always output local timers' values on Proc 0");
    plist->set ("writeGlobalStats", writeGlobalStats, "Always output global "
                "statistics, even if there is only one process in the "
                "communicator");
    plist->set ("writeZeroTimers", writeZeroTimers, "Generate output for "
                "timers that have never been called");
    return rcp_const_cast<const ParameterList> (plist);
  }

  void
  TimeMonitor::setReportParameters (const RCP<ParameterList>& params)
  {
    ETimeMonitorReportFormat reportFormat = REPORT_FORMAT_TABLE;
    ECounterSetOp setOp = Intersection;
    bool alwaysWriteLocal = false;
    bool writeGlobalStats = true;
    bool writeZeroTimers = true;

    if (params.is_null ()) {
      // If we've set parameters before, leave their current values.
      // Otherwise, set defaults (below).
      if (setParams_) {
        return;
      }
    }
    else { // params is nonnull.  Let's read it!
      params->validateParametersAndSetDefaults (*getValidReportParameters ());

      reportFormat = getIntegralValue<ETimeMonitorReportFormat> (*params, "Report format");
      setOp = getIntegralValue<ECounterSetOp> (*params, "How to merge timer sets");
      alwaysWriteLocal = params->get<bool> ("alwaysWriteLocal");
      writeGlobalStats = params->get<bool> ("writeGlobalStats");
      writeZeroTimers = params->get<bool> ("writeZeroTimers");
    }
    // Defer setting state until here, to ensure the strong exception
    // guarantee for this method (either it throws with no externally
    // visible state changes, or it returns normally).
    reportFormat_ = reportFormat;
    setOp_ = setOp;
    alwaysWriteLocal_ = alwaysWriteLocal;
    writeGlobalStats_ = writeGlobalStats;
    writeZeroTimers_ = writeZeroTimers;

    setParams_ = true; // Yay, we successfully set parameters!
  }

  void
  TimeMonitor::report (Ptr<const Comm<int> > comm,
                       std::ostream& out,
                       const RCP<ParameterList>& params)
  {
    setReportParameters (params);

    if (reportFormat_ == REPORT_FORMAT_YAML) {
      summarizeToYaml (comm, out);
    }
    else if (reportFormat_ == REPORT_FORMAT_TABLE) {
      summarize (comm, out, alwaysWriteLocal_, writeGlobalStats_,
                 writeZeroTimers_, setOp_);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "TimeMonitor::report: "
        "Invalid report format.  This should never happen; ParameterList "
        "validation should have caught this.  Please report this bug to the "
        "Teuchos developers.");
    }
  }

  void
  TimeMonitor::report (std::ostream& out,
                       const RCP<ParameterList>& params)
  {
    RCP<const Comm<int> > comm = getDefaultComm ();
    report (comm.ptr (), out, params);
  }

} // namespace Teuchos
