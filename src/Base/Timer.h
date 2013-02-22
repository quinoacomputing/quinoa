//******************************************************************************
/*!
  \file      src/Base/Timer.h
  \author    J. Bakosi
  \date      Fri Feb 22 15:49:33 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Timer
  \details   Timer
*/
//******************************************************************************
#ifndef Timer_h
#define Timer_h

#include <string>
#include <chrono>

#include <QuinoaTypes.h>

namespace Quinoa {

using TimerIndex = int;

const TimerIndex MAX_TIMERS = 32;

using namespace std;

// Shorthands for durations
using hours = chrono::hours;
using minutes = chrono::minutes;
using seconds = chrono::seconds;

//! Time representation in hours : minutess : seconds
struct HMS {
  hours h;
  minutes m;
  seconds s;
};

//! Quinoa::Timer
class Timer {

  // Shorthand for clock
  using clock = chrono::high_resolution_clock;

  // Shorthand for float seconds
  using dsec = chrono::duration<real>;

  //! Timer struct
  struct Clock {
    string name;
    bool used;
    clock::time_point start;
    clock::time_point now;
  };

  public:
    //! Constructor
    Timer();

    //! Destructor
    virtual ~Timer() {}

    //! Create new timer
    TimerIndex create(const string& label);

    //! Start timer
    void start(const TimerIndex id);

    //! Return time elapsed between start and stop
    real query(const TimerIndex id);

    //! Return time elapsed between start and stop
    void query(const TimerIndex id, HMS& hms);

    //! Estimate time for accomplishment
    void eta(const TimerIndex id,
             const real term,
             const real time,
             const real dt,
             const int nstep,
             const int it,
             HMS& hms);

  private:
    //! Don't permit copy constructor
    Timer(const Timer&) = delete;
    //! Don't permit copy assigment
    Timer& operator=(const Timer&) = delete;
    //! Don't permit move constructor
    Timer(Timer&&) = delete;
    //! Don't permit move assigment
    Timer& operator=(Timer&&) = delete;

    Clock m_timer[MAX_TIMERS];      //!< Timers
};

} // namespace Quinoa

#endif // Timer_h
