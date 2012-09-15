//******************************************************************************
/*!
  \file      src/Base/Exception.h
  \author    J. Bakosi
  \date      Fri Sep 14 17:32:08 2012
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

    //! Copy operator
    Exception(const Exception&);

    //! Destructor
    ~Exception() {}

    //! Handle Exception
    ErrCode handleException(Driver* driver);

  private:
    //! Don't permit assigment operator
    Exception& operator=(const Exception&);

    //! Exception type (CUMULATIVE, WARNING, ERROR, etc.)
    ExceptType m_except;
};

} // namespace Quinoa

#endif // Exception_h
