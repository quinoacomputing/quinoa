//******************************************************************************
/*!
  \file      src/Base/TimerException.h
  \author    J. Bakosi
  \date      Fri Apr 26 15:01:45 2013
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
                            const unsigned int& line) :
      Exception(except, file, func, line), m_except(timerExcept) {}

    //! Move constructor for throws, default compiler generated
    TimerException(TimerException&&) = default;

    //! Destructor
    virtual ~TimerException() noexcept = default;

    //! Handle TimerException
    virtual ErrCode handleException(Driver* const driver);

  private:
    //! Don't permit copy constructor
    TimerException(const TimerException&) = delete;
    //! Don't permit copy assignment
    TimerException& operator=(const TimerException&) = delete;
    //! Don't permit move assignment
    TimerException& operator=(TimerException&&) = delete;

    //! Timer exception type (BAD_ALLOC, BAD_INSERT, BAD_NAME, etc.)
    const TimerExceptType m_except;
};

} // namespace Quinoa

#endif // TimerException_h
