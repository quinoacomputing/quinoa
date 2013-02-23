//******************************************************************************
/*!
  \file      src/Base/Timer.C
  \author    J. Bakosi
  \date      Sat 23 Feb 2013 08:25:31 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Timer
  \details   Timer
*/
//******************************************************************************

#include <iostream>

#include <Timer.h>
#include <TimerException.h>

using namespace Quinoa;

Timer::Timer()
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  Clock zero = { "", false, clock::now(), clock::now() };

  for (int i=0; i<MAX_TIMERS; ++i) {
    m_timer[i] = zero;
  }
}

TimerIdx
Timer::create(const string& label)
//******************************************************************************
//  Create new timer
//! \param[in]  label  Name of timer
//! \return            Timer index
//! \author J. Bakosi
//******************************************************************************
{
  Assert(label.size() > 0, TimerException,FATAL,EMPTY_CLOCK_NAME);

  // Find an unused timer
  bool found = false;
  TimerIdx id=0, i=0;
  while (!found && i<MAX_TIMERS) {
    if (!m_timer[i].used) {
      found = true;
      id = i;
    }
    ++i;
  }

  Assert(found, TimerException,FATAL,TOO_MANY_TIMERS);

  m_timer[id].name = label;
  m_timer[id].used = true;

  return id;
}

void
Timer::start(const TimerIdx id)
//******************************************************************************
//  Start timer
//! \param[in]  id     Timer index
//! \author J. Bakosi
//******************************************************************************
{
  m_timer[id].start = clock::now();
}

real
Timer::query(const TimerIdx id)
//******************************************************************************
//  Return time elapsed between start and stop for timer as real
//! \param[in]  id     Timer index
//! \author J. Bakosi
//******************************************************************************
{
  using namespace std::chrono;

  // Get time stamp
  m_timer[id].now = clock::now();

  // Compute time difference between start and now
  dsec elapsed = duration_cast<dsec>(m_timer[id].now - m_timer[id].start);

  // Return elapsed time
  return elapsed.count();
}

void
Timer::query(const TimerIdx id, Watch& watch)
//******************************************************************************
//  Return time elapsed between start and stop for timer as h:m:s
//! \param[in]  id     Timer index
//! \param[out] watch  Watch holding time in hours:minutes:seconds
//! \author J. Bakosi
//******************************************************************************
{
  using namespace std::chrono;

  // Get time stamp
  m_timer[id].now = clock::now();

  // Compute time difference between start and now
  dsec elapsed = m_timer[id].now - m_timer[id].start;

  // Put elapsed time in watch as hours:minutes:seconds
  watch.h = duration_cast<hours>(elapsed);
  watch.m = duration_cast<minutes>(elapsed);
  watch.s = duration_cast<seconds>(elapsed);
}

void
Timer::eta(const TimerIdx id,
           const real term,
           const real time,
           const int nstep,
           const int it,
           Watch& elapsedWatch,
           Watch& estimatedWatch)
//******************************************************************************
//  Estimate time for accomplishment
//! \param[in]  id              Timer index
//! \param[in]  term            Time to terminate time stepping
//! \param[in]  time            Current time
//! \param[in]  nstep           Number time steps to take
//! \param[in]  it              Current iteration
//! \param[out] elapsedWatch    Elapsed time in h:m:s
//! \param[out] etimatedWatch   Estimated time for accomplishmet in h:m:s
//! \author J. Bakosi
//******************************************************************************
{
  using namespace std::chrono;

  dsec elapsed, estimated;

  if (it == 0) {
    // First iteration, just return zero
    elapsed = estimated = clock::duration::zero();
  } else {
    // Get time stamp
    m_timer[id].now = clock::now();

    // Compute time difference between start and now
    elapsed = m_timer[id].now - m_timer[id].start;

    // Estimate time until term in seconds
    dsec est_term = elapsed * (term-time) / time;
    // Estimate time until nstep in seconds
    dsec est_nstep = elapsed * static_cast<real>(nstep-it) / it;

    // Time stepping will stop at term or nstep, whichever is sooner
    estimated = min(est_term, est_nstep);
  }

  // Put elapsed time in watch as hours:minutes:seconds
  estimatedWatch.h = duration_cast<hours>(estimated);
  estimatedWatch.m = duration_cast<minutes>(estimated);
  estimatedWatch.s = duration_cast<seconds>(estimated);
  // Put estimated time in watch as hours:minutes:seconds
  elapsedWatch.h = duration_cast<hours>(elapsed);
  elapsedWatch.m = duration_cast<minutes>(elapsed);
  elapsedWatch.s = duration_cast<seconds>(elapsed);
}
