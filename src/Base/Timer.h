//******************************************************************************
/*!
  \file      src/Base/Timer.h
  \author    J. Bakosi
  \date      Thu 21 Feb 2013 10:31:48 PM MST
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

//! Quinoa::Timer
class Timer {

  using clock = std::chrono::high_resolution_clock;

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
