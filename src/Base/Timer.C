// *****************************************************************************
/*!
  \file      src/Base/Timer.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Timer definition
  \details   Timer definition. Timer is a simple class to do timing various
    parts of the code in a portable way. The functionality is intended to be
    very minimal and simple, but still convenient to use, with as little state
    as possible. For an example client code, see Main.
*/
// *****************************************************************************

#include <iostream>

#include "Timer.h"

using tk::Timer;

tk::Timer::Watch
Timer::hms() const
// *****************************************************************************
//  Return time elapsed between start and stop for timer as hours, minutes, and
//  seconds.
//! \return Time elapsed between start and stop as hours, minutes, and seconds,
//!   as a Watch struct.
//! \author J. Bakosi
// *****************************************************************************
{
  using std::chrono::duration_cast;

  // Compute time difference between start and now in seconds
  Dsec elapsed = clock::now() - m_start;

  // Put elapsed time in watch as hours:minutes:seconds
  Watch watch( duration_cast< hours >( elapsed ),
               duration_cast< minutes >( elapsed ) % hours(1),
               duration_cast< seconds >( elapsed ) % minutes(1) );
  return watch;
}

void
Timer::eta( tk::real term, tk::real time, uint64_t nstep, uint64_t it,
            Watch& elapsedWatch, Watch& estimatedWatch ) const
// *****************************************************************************
//  Estimate time for accomplishment
//! \param[in]  term            Time at which to terminate time stepping
//! \param[in]  time            Current time
//! \param[in]  nstep           Max number of time steps to take
//! \param[in]  it              Current iteration count
//! \param[out] elapsedWatch    Elapsed time in h:m:s
//! \param[out] estimatedWatch  Estimated time for accomplishmet in h:m:s
//! \author J. Bakosi
// *****************************************************************************
{
  using std::chrono::duration_cast;

  Dsec elapsed, estimated;

  if (it == 0) {

    // First iteration, just return zero
    elapsed = estimated = clock::duration::zero();

  } else {

    // Compute time difference between start and now in seconds
    elapsed = clock::now() - m_start;

    // Estimate time until term in seconds
    Dsec est_term = elapsed * (term-time) / time;
    // Estimate time until nstep in seconds
    Dsec est_nstep = elapsed * static_cast<tk::real>(nstep-it) / it;

    // Time stepping will stop at term or nstep, whichever is sooner
    estimated = min(est_term, est_nstep);

  }

  // Put elapsed time in watch as hours:minutes:seconds
  elapsedWatch.hrs = duration_cast< hours >( elapsed );
  elapsedWatch.min = duration_cast< minutes >( elapsed ) % hours(1);
  elapsedWatch.sec = duration_cast< seconds >( elapsed ) % minutes(1);
  // Put estimated time in watch as hours:minutes:seconds
  estimatedWatch.hrs = duration_cast< hours >( estimated );
  estimatedWatch.min = duration_cast< minutes >( estimated ) % hours(1);
  estimatedWatch.sec = duration_cast< seconds >( estimated ) % minutes(1);
}

namespace tk {

Timer::Watch
hms( tk::real stamp )
// *****************************************************************************
//! Convert existing time stamp as a real to Watch (global-scope)
//! \param[in] stamp Time stamp as a real number
//! \return Time as hours, minutes, and seconds, as a Watch struct.
//! \author J. Bakosi
// *****************************************************************************
{
  using std::chrono::duration_cast;
  const auto d = Timer::Dsec( stamp );
  return
    Timer::Watch( duration_cast< Timer::hours >( d ),
                  duration_cast< Timer::minutes >( d ) % Timer::hours(1),
                  duration_cast< Timer::seconds >( d ) % Timer::minutes(1) );
}

} // tk::
