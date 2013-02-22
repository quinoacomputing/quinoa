//******************************************************************************
/*!
  \file      src/Base/TimerException.h
  \author    J. Bakosi
  \date      Thu 21 Feb 2013 09:06:43 PM MST
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
    TimerException(ExceptType except,
                   TimerExceptType timerExcept,
                   const string& file,
                   const string& func,
                   const unsigned int& line) :
      Exception(except, file, func, line), m_except(timerExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    TimerException(TimerException&&) = default;

    //! Destructor
    //virtual ~TimerException() {}

    //! Handle TimerException
    virtual ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy constructor
    TimerException(const TimerException&) = delete;
    //! Don't permit copy assignment
    TimerException& operator=(const TimerException&) = delete;
    //! Don't permit move assignment
    TimerException& operator=(TimerException&&) = delete;

    //! Timer exception type (BAD_ALLOC, BAD_INSERT, BAD_NAME, etc.)
    TimerExceptType m_except;
};

} // namespace Quinoa

#endif // TimerException_h
