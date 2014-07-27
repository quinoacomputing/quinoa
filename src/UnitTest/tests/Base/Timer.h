//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/Timer.h
  \author    J. Bakosi
  \date      Sat 26 Jul 2014 09:33:27 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Unit tests for Base/Timer.h
  \details   Unit tests for Base/Timer.h
*/
//******************************************************************************
#ifndef test_Timer_h
#define test_Timer_h

#include <unistd.h>
#include <tut/tut.hpp>
#include <Timer.h>

namespace tut {

//! All tests in group inherited from this base
struct Timer_common {};

//! Test group shortcuts
using Timer_group = test_group< Timer_common >;
using Timer_object = Timer_group::object;

//! Define test group
Timer_group Timer( "Base/Timer" );

//! Test definitions for group

//! Test timing a 0.1s duration as float with 1/10th millisecond precision
template<> template<>
void Timer_object::test< 1 >() {
  set_test_name( "measure 0.1s using dsec() with 1/10th ms prec" );

  tk::Timer timer;
  usleep( 100000 );    // in micro-seconds, sleep for 0.1 second
  // test if time measured with at least 1/10th of a millisecond precision
  ensure_equals( "time 0.1s elapsed as float", timer.dsec(), 0.1, 1.0e-4 );
}

//! Test timing a 0.1 duration as h:m:s with 1/10th millisecond precision
template<> template<>
void Timer_object::test< 2 >() {
  set_test_name( "measure 0.1s using hms() with 1/10th ms prec" );

  tk::Timer timer;
  usleep( 100000 );    // in micro-seconds, sleep for 0.1 second
  const auto stamp = timer.hms();
  // test if time measured with at least 1/10th of a millisecond precision
  ensure_equals( "time 0.1s elapsed as hrs", stamp.hrs.count(), 0.0, 1.0e-4 );
  ensure_equals( "time 0.1s elapsed as min", stamp.min.count(), 0.0, 1.0e-4 );
  ensure_equals( "time 0.1s elapsed as sec", stamp.sec.count(), 0.0, 1.0e-4 );
}

//! Test estimated time elapsed and to accomplishment triggered by term
template<> template<>
void Timer_object::test< 3 >() {
  set_test_name( "ETE and ETA triggered by terminate time" );

  // Setup a duration case
  tk::real term = 5.0;      // time at which to terminate time stepping
  tk::real time = 1.0;      // current time
  uint64_t nstep = 1000;    // max number of time steps to take (large)
  uint64_t it = 1;          // current iteration

  tk::Timer timer;
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  tk::Timer::Watch ete, eta;
  timer.eta( term, time, nstep, it, ete, eta );
  // test estimated time elapsed with 1/10th ms precision
  ensure_equals( "estimated time elapsed in hrs", ete.hrs.count(),
                 0.0, 1.0e-4 );
  ensure_equals( "estimated time elapsed in min", ete.min.count(),
                 0.0, 1.0e-4 );
  ensure_equals( "estimated time elapsed in sec", ete.sec.count(),
                 1.0, 1.0e-4 );
  // test estimated time to accomplishment with 1/10th ms precision
  ensure_equals( "estimated time to accomlishment in hrs", eta.hrs.count(),
                 0.0, 1.0e-4 );
  ensure_equals( "estimated time to accomlishment in min", eta.min.count(),
                 0.0, 1.0e-4 );
  ensure_equals( "estimated time to accomlishment in sec", eta.sec.count(),
                 4.0, 1.0e-4 );
}

//! Test estimated time elapsed and to accomplishment triggered by nstep
template<> template<>
void Timer_object::test< 4 >() {
  set_test_name( "ETE and ETA triggered by max number of steps" );

  // Setup a duration case
  tk::real term = 500.0;    // time at which to terminate time stepping (large)
  tk::real time = 1.0;      // current time
  uint64_t nstep = 100;     // max number of time steps to take
  uint64_t it = 1;          // current iteration

  tk::Timer timer;
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  tk::Timer::Watch ete, eta;
  timer.eta( term, time, nstep, it, ete, eta );
  // test estimated time elapsed with 1/10th ms precision
  ensure_equals( "estimated time elapsed in hrs", ete.hrs.count(), 0.0, 1.0e-4 );
  ensure_equals( "estimated time elapsed in min", ete.min.count(), 0.0, 1.0e-4 );
  ensure_equals( "estimated time elapsed in sec", ete.sec.count(), 1.0, 1.0e-4 );
  // test estimated time to accomplishment with 1/10th ms precision
  ensure_equals( "estimated time to accomlishment in hrs", eta.hrs.count(),
                 0.0, 1.0e-4 );
  ensure_equals( "estimated time to accomlishment in min", eta.min.count(),
                 1.0, 1.0e-4 );
  ensure_equals( "estimated time to accomlishment in sec", eta.sec.count(),
                 39.0, 1.0e-4 );
}

} // tut::

#endif // test_Timer_h
