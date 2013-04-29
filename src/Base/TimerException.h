//******************************************************************************
/*!
  \file      src/Base/TimerException.h
  \author    J. Bakosi
  \date      Mon Apr 29 15:57:15 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TimerException class declaration
  \details   TimerException class declaration
*/
//******************************************************************************
#ifndef TimerException_h
#define TimerException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! Timer exception types
enum TimerExceptType { EMPTY_CLOCK_NAME=0,     //!< Noname timer
                       TOO_MANY_TIMERS,        //!< Too many timers registered
                       NUM_TIMER_EXCEPT
};

//! Timer exception error messages
const string TimerMsg[NUM_TIMER_EXCEPT] = {
  "Must specify a name with non-zero length",
  "Too many timers registered"
};

//! TimerException : Exception
class TimerException : public Exception {

  public:
    //! Constructor
    explicit TimerException(const ExceptType except,
                            const TimerExceptType timerExcept,
                            const string& file,
                            const string& func,
                            unsigned int line) noexcept :
      Exception(except,
                TimerMsg[static_cast<int>(timerExcept)],
                file,
                func,
                line) {}

    //! Move constructor for throws, default compiler generated
    TimerException(TimerException&&) = default;

    //! Destructor
    virtual ~TimerException() noexcept = default;

  private:
    //! Don't permit copy constructor
    TimerException(const TimerException&) = delete;
    //! Don't permit copy assignment
    TimerException& operator=(const TimerException&) = delete;
    //! Don't permit move assignment
    TimerException& operator=(TimerException&&) = delete;
};

} // namespace Quinoa

#endif // TimerException_h
