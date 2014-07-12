//******************************************************************************
/*!
  \file      src/Base/Timer.h
  \author    J. Bakosi
  \date      Tue 01 Jul 2014 08:38:20 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Timer
  \details   Timer
*/
//******************************************************************************
#ifndef Timer_h
#define Timer_h

#include <chrono>

#include <Types.h>

namespace tk {

//! Watch stores time in hours:minutes:seconds
struct Watch {
  std::chrono::hours hrs;
  std::chrono::minutes min;
  std::chrono::seconds sec;
  Watch( std::chrono::hours&& h,
         std::chrono::minutes&& m,
         std::chrono::seconds&& s ) :
    hrs( std::move(h) ), min( std::move(m) ), sec( std::move(s) ) {}
};

//! Timer
class Timer {

  private:
    // Shorthand for clock, set clock type
    using clock = std::chrono::high_resolution_clock;

    // Shorthand for real duration
    using Dsec = std::chrono::duration< real >;

  public:
    //! Constructor
    Timer() : m_start( clock::now() ) {}

    //! Return time elapsed between start and stop as a real number
    real dsec() const {
      return std::chrono::duration_cast< Dsec >(clock::now() - m_start).count();
    }

    //! Return time elapsed between start and stop as h:m:s
    Watch hms() const;

    //! Estimate time for accomplishment
    void eta( const real term,
              const real time,
              const uint64_t nstep,
              const uint64_t it,
              Watch& elapsedWatch,
              Watch& estimatedWatch ) const;

  private:
    clock::time_point m_start;  //!< Time stamp at start
};

} // tk::

#endif // Timer_h
