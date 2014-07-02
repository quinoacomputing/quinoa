//******************************************************************************
/*!
  \file      src/Base/Timer.C
  \author    J. Bakosi
  \date      Tue 01 Jul 2014 08:38:32 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Timer
  \details   Timer
*/
//******************************************************************************

#include <iostream>

#include <Timer.h>

using tk::Timer;
using tk::Watch;

Watch
Timer::hms() const
//******************************************************************************
//  Return time elapsed between start and stop for timer as h:m:s
//! \param[out] watch  Watch holding time in hours:minutes:seconds
//! \author J. Bakosi
//******************************************************************************
{
  using std::chrono::duration_cast;
  using std::chrono::hours;
  using std::chrono::minutes;
  using std::chrono::seconds;

  // Compute time difference between start and now in seconds
  Dsec elapsed = (clock::now() - m_start) / 1000.0;

  // Put elapsed time in watch as hours:minutes:seconds
  Watch watch( duration_cast< hours >( elapsed ),
               duration_cast< minutes >( elapsed ) % hours(1),
               duration_cast< seconds >( elapsed ) % minutes(1) );
  return watch;
}

void
Timer::eta( const tk::real term,
            const tk::real time,
            const uint64_t nstep,
            const uint64_t it,
            Watch& elapsedWatch,
            Watch& estimatedWatch ) const
//******************************************************************************
//  Estimate time for accomplishment
//! \param[in]  term            Time to terminate time stepping
//! \param[in]  time            Current time
//! \param[in]  nstep           Number time steps to take
//! \param[in]  it              Current iteration
//! \param[out] elapsedWatch    Elapsed time in h:m:s
//! \param[out] estimatedWatch  Estimated time for accomplishmet in h:m:s
//! \author J. Bakosi
//******************************************************************************
{
  using std::chrono::duration_cast;
  using std::chrono::hours;
  using std::chrono::minutes;
  using std::chrono::seconds;

  Dsec elapsed, estimated;

  if (it == 0) {

    // First iteration, just return zero
    elapsed = estimated = clock::duration::zero();

  } else {

    // Compute time difference between start and now in seconds
    elapsed = (clock::now() - m_start) / 1000.0;

    // Estimate time until term in seconds
    Dsec est_term = elapsed * (term-time) / time;
    // Estimate time until nstep in seconds
    Dsec est_nstep = elapsed * static_cast<tk::real>(nstep-it) / it;

    // Time stepping will stop at term or nstep, whichever is sooner
    estimated = min(est_term, est_nstep);

  }

  // Put elapsed time in watch as hours:minutes:seconds
  elapsedWatch.hrs = duration_cast<hours>(elapsed);
  elapsedWatch.min = duration_cast<minutes>(elapsed) % hours(1);
  elapsedWatch.sec = duration_cast<seconds>(elapsed) % minutes(1);
  // Put estimated time in watch as hours:minutes:seconds
  estimatedWatch.hrs = duration_cast<hours>(estimated);
  estimatedWatch.min = duration_cast<minutes>(estimated) % hours(1);
  estimatedWatch.sec = duration_cast<seconds>(estimated) % minutes(1);
}
