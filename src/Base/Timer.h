//******************************************************************************
/*!
  \file      src/Base/Timer.h
  \author    J. Bakosi
  \date      Wed Apr 23 16:40:19 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Timer
  \details   Timer
*/
//******************************************************************************
#ifndef Timer_h
#define Timer_h

#include <string>
#include <chrono>
#include <vector>

#include <Types.h>

namespace tk {

using TimerId = std::size_t;

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
      clock::time_point start;  //!< Time stamp at start
    };

  public:
    //! Constructor
    explicit Timer() = default;

    //! Destructor
    virtual ~Timer() noexcept = default;

    //! Create new timer
    TimerId create( std::string label ) const;

    //! Create new timer
    TimerId create( const char* label ) const {
      return create( std::string(label) );
    }

    //! Start timer
    void start( TimerId id ) const;

    //! Start all timer
    void start() const;

    //! Return time elapsed between start and stop
    real query( TimerId id ) const;

    //! Return time elapsed between start and stop
    void query( TimerId id, Watch& watch ) const;

    //! Estimate time for accomplishment
    void eta( TimerId id,
              const real term,
              const real time,
              const uint64_t nstep,
              const uint64_t it,
              Watch& elapsedWatch,
              Watch& estimatedWatch ) const;

    //! Return timer name
    std::string name( TimerId id ) const;

  private:
    //! Don't permit copy constructor
    Timer(const Timer&) = delete;
    //! Don't permit copy assigment
    Timer& operator=(const Timer&) = delete;
    //! Don't permit move constructor
    Timer(Timer&&) = delete;
    //! Don't permit move assigment
    Timer& operator=(Timer&&) = delete;

    mutable std::vector< Clock > m_timer;      //!< Timers
};

} // tk::

#endif // Timer_h
