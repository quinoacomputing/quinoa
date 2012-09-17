//******************************************************************************
/*!
  \file      src/Base/Exception.h
  \author    J. Bakosi
  \date      Sun 16 Sep 2012 05:54:34 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class declaration
  \details   Exception base class declaration
*/
//******************************************************************************
#ifndef Exception_h
#define Exception_h

#include <QuinoaTypes.h>

namespace Quinoa {

//! Exception types
enum class ExceptType { CUMULATIVE=0,  //!< Only several will produce a warning
                        WARNING,       //!< Warning: output message
                        ERROR,         //!< Error: output but will not interrupt
                        UNCAUGHT,      //!< Uncaught: will interrupt
                        FATAL,         //!< Fatal error: will interrupt
                        NUM_EXCEPT
};
//! Number of exception types
const Int NUM_EXCEPT = static_cast<Int>(ExceptType::NUM_EXCEPT);

//! Error codes for the OS (or whatever calls Quinoa)
enum class ErrCode { NO_ERROR=0,       //!< Everything went fine
                     NONFATAL,         //!< Exception occurred but continue
                     FATAL,            //!< Fatal error occurred
                     NUM_ERR_CODE
};
//! Number of error codes
const Int NUM_ERR_CODE = static_cast<Int>(ErrCode::NUM_ERR_CODE);

class Driver;

//! Exception base class
class Exception {

  public:
    //! Constructor
    Exception(ExceptType except) : m_except(except) {}

    //! Handle Exception
    ErrCode handleException(Driver* driver);

  protected:
    //! Move constructor, necessary for throws, default compiler generated,
    //! can only be thrown from within derived Exception classes
    Exception(Exception&&) = default;

    //! Destructor
    ~Exception() = default;

  private:
    //! Don't permit copy constructor
    Exception(const Exception&) = delete;
    //! Don't permit copy assignment
    Exception& operator=(const Exception&) = delete;
    //! Don't permit move assignment
    Exception& operator=(Exception&&) = delete;

    //! Exception type (CUMULATIVE, WARNING, ERROR, etc.)
    ExceptType m_except;
};

} // namespace Quinoa

#endif // Exception_h
