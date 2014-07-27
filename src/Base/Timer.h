//******************************************************************************
/*!
  \file      src/Base/Timer.h
  \author    J. Bakosi
  \date      Sat 26 Jul 2014 09:39:31 PM MDT
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

//! Timer
class Timer {

  private:
    using Dsec = std::chrono::duration< real >;
    using hours = std::chrono::hours;
    using minutes = std::chrono::minutes;
    using seconds = std::chrono::seconds;

  public:
    // Shorthand for clock, setting clock type
    using clock = std::chrono::high_resolution_clock;

    //! Watch stores time in hours:minutes:seconds
    struct Watch {
      hours hrs;
      minutes min;
      seconds sec;
      //! Zero constructor
      explicit Watch() :
        hrs( std::chrono::duration_cast< hours >( clock::duration::zero()) ),
        min( std::chrono::duration_cast< minutes >( clock::duration::zero() ) ),
        sec( std::chrono::duration_cast< seconds >( clock::duration::zero() ) )
      {}
      //! Fill constructor
      explicit Watch( hours&& h, minutes&& m, seconds&& s ) :
        hrs( std::move(h) ), min( std::move(m) ), sec( std::move(s) ) {}
    };

    //! Constructor
    explicit Timer() : m_start( clock::now() ) {}

    //! Return time elapsed between start and stop as a real number
    real dsec() const {
      return std::chrono::duration_cast< Dsec >(clock::now() - m_start).count();
    }

    //! Return time elapsed between start and stop as h:m:s
    Watch hms() const;

    //! Estimate time for accomplishment
    void eta( real term, real time, uint64_t nstep, uint64_t it,
              Watch& elapsedWatch, Watch& estimatedWatch ) const;

  private:
    clock::time_point m_start;  //!< Time stamp at start
};

} // tk::

#endif // Timer_h
