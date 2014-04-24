//******************************************************************************
/*!
  \file      src/Base/Timer.C
  \author    J. Bakosi
  \date      Thu 24 Apr 2014 06:49:37 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Timer
  \details   Timer
*/
//******************************************************************************

#include <iostream>

#include <Timer.h>
#include <Exception.h>

using tk::Timer;

tk::TimerId
Timer::create( std::string label ) const
//******************************************************************************
//  Create new timer
//! \param[in]  label  Name of timer
//! \return            Timer index
//! \author J. Bakosi
//******************************************************************************
{
  Assert( label.size() > 0, ExceptType::FATAL,
          "Must give a non-empty string as timer name" );

  m_timer.push_back( { std::move(label), clock::now() } );

  return m_timer.size()-1;
}

void
Timer::start( TimerId id ) const
//******************************************************************************
//  Start timer
//! \param[in]  id     Timer index
//! \author J. Bakosi
//******************************************************************************
{
  Assert( m_timer.size() > id, ExceptType::FATAL,
          "Attempt to start a non-existent timer" );

  m_timer[id].start = clock::now();
}

void
Timer::start() const
//******************************************************************************
//  Start all registered timers
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !m_timer.empty(), ExceptType::WARNING,
          "There are not registered timers" );

  for (auto& t : m_timer) t.start = clock::now();
}

tk::real
Timer::query( TimerId id ) const
//******************************************************************************
//  Return time elapsed between start and stop for timer as real
//! \param[in]  id     Timer index
//! \author J. Bakosi
//******************************************************************************
{
  Assert( m_timer.size() > id, ExceptType::FATAL,
          "Attempt to query a non-existent timer" );

  // Compute time difference between start and now in real
  dsec elapsed =
    std::chrono::duration_cast< dsec >( clock::now() - m_timer[id].start );

  // Return elapsed time
  return elapsed.count();
}

void
Timer::query( TimerId id, Watch& watch ) const
//******************************************************************************
//  Return time elapsed between start and stop for timer as h:m:s
//! \param[in]  id     Timer index
//! \param[out] watch  Watch holding time in hours:minutes:seconds
//! \author J. Bakosi
//******************************************************************************
{
  Assert( m_timer.size() > id, ExceptType::FATAL,
          "Attempt to query a non-existent timer" );

  using std::chrono::duration_cast;
  using std::chrono::hours;
  using std::chrono::minutes;
  using std::chrono::seconds;

  // Compute time difference between start and now in seconds
  dsec elapsed = (clock::now() - m_timer[id].start) / 1000.0;

  // Put elapsed time in watch as hours:minutes:seconds
  watch.h = duration_cast< hours >( elapsed );
  watch.m = duration_cast< minutes >( elapsed ) % hours(1);
  watch.s = duration_cast< seconds >( elapsed ) % minutes(1);
}

void
Timer::eta( TimerId id,
            const tk::real term,
            const tk::real time,
            const uint64_t nstep,
            const uint64_t it,
            Watch& elapsedWatch,
            Watch& estimatedWatch ) const
//******************************************************************************
//  Estimate time for accomplishment
//! \param[in]  id              Timer index
//! \param[in]  term            Time to terminate time stepping
//! \param[in]  time            Current time
//! \param[in]  nstep           Number time steps to take
//! \param[in]  it              Current iteration
//! \param[out] elapsedWatch    Elapsed time in h:m:s
//! \param[out] estimatedWatch  Estimated time for accomplishmet in h:m:s
//! \author J. Bakosi
//******************************************************************************
{
  Assert( m_timer.size() > id, ExceptType::FATAL,
          "Attempt to query a non-existent timer" );

  using std::chrono::duration_cast;
  using std::chrono::hours;
  using std::chrono::minutes;
  using std::chrono::seconds;

  dsec elapsed, estimated;

  if (it == 0) {

    // First iteration, just return zero
    elapsed = estimated = clock::duration::zero();

  } else {

    // Compute time difference between start and now in seconds
    elapsed = (clock::now() - m_timer[id].start) / 1000.0;

    // Estimate time until term in seconds
    dsec est_term = elapsed * (term-time) / time;
    // Estimate time until nstep in seconds
    dsec est_nstep = elapsed * static_cast<tk::real>(nstep-it) / it;

    // Time stepping will stop at term or nstep, whichever is sooner
    estimated = min(est_term, est_nstep);

  }

  // Put elapsed time in watch as hours:minutes:seconds
  elapsedWatch.h = duration_cast<hours>(elapsed);
  elapsedWatch.m = duration_cast<minutes>(elapsed) % hours(1);
  elapsedWatch.s = duration_cast<seconds>(elapsed) % minutes(1);
  // Put estimated time in watch as hours:minutes:seconds
  estimatedWatch.h = duration_cast<hours>(estimated);
  estimatedWatch.m = duration_cast<minutes>(estimated) % hours(1);
  estimatedWatch.s = duration_cast<seconds>(estimated) % minutes(1);
}

std::string
Timer::name( TimerId id ) const
//******************************************************************************
//  Return timer name
//! \param[in]  id     Timer index
//! \author J. Bakosi
//******************************************************************************
{
  Assert( m_timer.size() > id, ExceptType::FATAL,
          "Attempt to query a non-existent timer" );

  return m_timer[id].name;
}
