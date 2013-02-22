//******************************************************************************
/*!
  \file      src/Base/Timer.C
  \author    J. Bakosi
  \date      Thu 21 Feb 2013 10:31:53 PM MST
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
//  Return time elapsed between start and stop for timer
//! \param[in]  id     Timer ID
//! \author J. Bakosi
//******************************************************************************
{
  using namespace std::chrono;

  m_timer[id].now = clock::now();

  duration<real> elapsed =
    duration_cast<duration<real>>(m_timer[id].now - m_timer[id].start);

  return elapsed.count();
}
