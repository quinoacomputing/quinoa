//******************************************************************************
/*!
  \file      src/Base/Timer.h
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 05:45:56 PM MDT
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

namespace quinoa {

using TimerIdx = int;

const TimerIdx MAX_TIMERS = 32;

//! Watch stores time in hours:minutes:seconds
struct Watch {
  std::chrono::hours h;
  std::chrono::minutes m;
  std::chrono::seconds s;
};

//! Timer
class Timer {

  private:
    // Shorthand for clock, set clock type
    using clock = std::chrono::high_resolution_clock;

    // Shorthand for float seconds
    using dsec = std::chrono::duration<real>;

    //! Timer struct
    struct Clock {
      std::string name;         //!< Timer name
      bool used;                //!< In use or not
      clock::time_point start;  //!< Time stamp at start
      clock::time_point now;    //!< Time stamp at a later time
    };

  public:
    //! Constructor
    explicit Timer();

    //! Destructor
    virtual ~Timer() noexcept = default;

    //! Create new timer
    TimerIdx create(const std::string& label) const;

    //! Start timer
    void start(const TimerIdx id) const;

    //! Return time elapsed between start and stop
    real query(const TimerIdx id) const;

    //! Return time elapsed between start and stop
    void query(const TimerIdx id, Watch& watch) const;

    //! Estimate time for accomplishment
    void eta(const TimerIdx id,
             const real term,
             const real time,
             const uint64_t nstep,
             const uint64_t it,
             Watch& elapsedWatch,
             Watch& estimatedWatch) const;

  private:
    //! Don't permit copy constructor
    Timer(const Timer&) = delete;
    //! Don't permit copy assigment
    Timer& operator=(const Timer&) = delete;
    //! Don't permit move constructor
    Timer(Timer&&) = delete;
    //! Don't permit move assigment
    Timer& operator=(Timer&&) = delete;

    mutable Clock m_timer[MAX_TIMERS];      //!< Timers
};

} // namespace quinoa

#endif // Timer_h
