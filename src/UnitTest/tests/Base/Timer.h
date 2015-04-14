//******************************************************************************
/*!
  \file      src/UnitTest/tests/Base/Timer.h
  \author    J. Bakosi
  \date      Tue 14 Apr 2015 11:45:44 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
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
struct Timer_common {
  double precision = 1.0e-3;    // required precision in seconds for timings
};

//! Test group shortcuts
using Timer_group = test_group< Timer_common, MAX_TESTS_IN_GROUP >;
using Timer_object = Timer_group::object;

//! Define test group
Timer_group Timer( "Base/Timer" );

//! Test definitions for group

//! Test timing a 0.1s duration as float with given precision
template<> template<>
void Timer_object::test< 1 >() {
  set_test_name( "measure 0.1s using dsec() with " + std::to_string(precision) +
                 "s prec" );

  tk::Timer timer;
  usleep( 100000 );    // in micro-seconds, sleep for 0.1 second
  // test if time measured with at least 1/10th of a millisecond precision
  ensure_equals( "time 0.1s elapsed as float", timer.dsec(), 0.1, precision );
}

//! Test timing a 1.0 duration as h:m:s with given precision
template<> template<>
void Timer_object::test< 2 >() {
  set_test_name( "measure 1.0s using hms() with " + std::to_string(precision) +
                 "s prec" );

  tk::Timer timer;
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  const auto stamp = timer.hms();
  // test if time measured with at least 1/10th of a millisecond precision
  ensure_equals( "time 1.0s elapsed as hrs",
                 static_cast<tk::real>(stamp.hrs.count()), 0.0, precision );
  ensure_equals( "time 1.0s elapsed as min",
                 static_cast<tk::real>(stamp.min.count()), 0.0, precision );
  ensure_equals( "time 1.0s elapsed as sec",
                 static_cast<tk::real>(stamp.sec.count()), 1.0, precision );
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
  // test estimated time elapsed with given precision
  ensure_equals( "estimated time elapsed in hrs",
                 static_cast<tk::real>(ete.hrs.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in min",
                 static_cast<tk::real>(ete.min.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in sec",
                 static_cast<tk::real>(ete.sec.count()), 1.0, precision );
  // test estimated time to accomplishment with given precision
  ensure_equals( "estimated time to accomlishment in hrs",
                 static_cast<tk::real>(eta.hrs.count()), 0.0, precision );
  ensure_equals( "estimated time to accomlishment in min",
                 static_cast<tk::real>(eta.min.count()), 0.0, precision );
  ensure_equals( "estimated time to accomlishment in sec",
                 static_cast<tk::real>(eta.sec.count()), 4.0, precision );
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
  // test estimated time elapsed with given precision
  ensure_equals( "estimated time elapsed in hrs",
                 static_cast<tk::real>(ete.hrs.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in min",
                 static_cast<tk::real>(ete.min.count()), 0.0, precision );
  ensure_equals( "estimated time elapsed in sec",
                 static_cast<tk::real>(ete.sec.count()), 1.0, precision );
  // test estimated time to accomplishment with given precision
  ensure_equals( "estimated time to accomlishment in hrs",
                 static_cast<tk::real>(eta.hrs.count()), 0.0, precision );
  ensure_equals( "estimated time to accomlishment in min",
                 static_cast<tk::real>(eta.min.count()), 1.0, precision );
  ensure_equals( "estimated time to accomlishment in sec",
                 static_cast<tk::real>(eta.sec.count()), 39.0, precision );
}

//! Test converting a 1.0s duration timed as a float to Timer::Watch
template<> template<>
void Timer_object::test< 5 >() {
  set_test_name( "convert time stamp in float to Watch" );

  tk::Timer timer;
  usleep( 1000000 );    // in micro-seconds, sleep for 1.0 second
  // convert time stamp in float to Timer::Watch
  const auto w = tk::hms( timer.dsec() );
  // test if time measured with at least 1/10th of a millisecond precision
  ensure_equals( "time 1.0s elapsed as float represented as Timer::Watch in hrs",
                 static_cast<tk::real>(w.hrs.count()), 0.0, precision );
  ensure_equals( "time 1.0s elapsed as float represented as Timer::Watch in min",
                 static_cast<tk::real>(w.min.count()), 0.0, precision );
  ensure_equals( "time 1.0s elapsed as float represented as Timer::Watch in sec",
                 static_cast<tk::real>(w.sec.count()), 1.0, precision );
}

} // tut::

#endif // test_Timer_h
