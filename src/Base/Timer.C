//******************************************************************************
/*!
  \file      src/Base/Timer.C
  \author    J. Bakosi
  \date      Fri Feb 22 16:52:14 2013
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

  for (int i=0; i<MAX_TIMERS; ++i)
    m_timer[i] = zero;
}

TimerIndex
Timer::create(const string& label)
//******************************************************************************
//  Create new timer
//! \param[in]  label  Name of timer
//! \return            Timer ID
//! \author J. Bakosi
//******************************************************************************
{
  Assert(label.size() > 0, TimerException,FATAL,EMPTY_CLOCK_NAME);

  // Find an unused timer
  bool found = false;
  TimerIndex id=0, i=0;
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
Timer::start(const TimerIndex id)
//******************************************************************************
//  Start timer
//! \param[in]  id     Timer ID
//! \author J. Bakosi
//******************************************************************************
{
  m_timer[id].start = clock::now();
}

real
Timer::query(const TimerIndex id)
//******************************************************************************
//  Return time elapsed between start and stop for timer as real
//! \param[in]  id     Timer ID
//! \author J. Bakosi
//******************************************************************************
{
  using namespace std::chrono;

  // Get time stamp
  m_timer[id].now = clock::now();

  // Compute time difference between start and now
  dsec elapsed = duration_cast<dsec>(m_timer[id].now - m_timer[id].start);

  // Return elapsed time in real
  return elapsed.count();
}

void
Timer::query(const TimerIndex id, HMS& hms)
//******************************************************************************
//  Return time elapsed between start and stop for timer as h:m:s
//! \param[in]  id     Timer ID
//! \param[out] hms    Struct HMS holding chrono::hours:minutes:seconds
//! \author J. Bakosi
//******************************************************************************
{
  using namespace std::chrono;

  // Get time stamp
  m_timer[id].now = clock::now();

  // Compute time difference between start and now
  dsec elapsed = m_timer[id].now - m_timer[id].start;

  // Put in elapsed time in hours : minutes : seconds
  hms = { duration_cast<hours>(elapsed),
          duration_cast<minutes>(elapsed),
          duration_cast<seconds>(elapsed) };
}

void
Timer::eta(const TimerIndex id,
           const real term,
           const real time,
           const real dt,
           const int nstep,
           const int it,
           HMS& hms)
//******************************************************************************
//  Estimate time for accomplishment
//! \param[in]  id     Timer ID
//! \author J. Bakosi
//******************************************************************************
{
  using namespace std::chrono;

  dsec estimated;

  if (it == 0) {
    // Will not divide by zero
    estimated = clock::duration::zero();
  } else {
    // Get time stamp
    m_timer[id].now = clock::now();

    // Compute time difference between start and now
    dsec elapsed = m_timer[id].now - m_timer[id].start;

    // Estimate time until term in seconds
    dsec est_term = elapsed * (term-time) / time;
    // Estimate time until nstep in seconds
    dsec est_nstep = elapsed * static_cast<real>(nstep-it) / it;

    // Time stepping will stop at whichever is sooner
    estimated = min(est_term, est_nstep);
  }

  // Put in estimated time in hours : minutes : seconds
  hms = { duration_cast<hours>(estimated),
          duration_cast<minutes>(estimated),
          duration_cast<seconds>(estimated) };
}
