//******************************************************************************
/*!
  \file      src/Base/Exception.h
  \author    J. Bakosi
  \date      Fri Sep 14 15:41:09 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class declaration
  \details   Exception base class declaration
*/
//******************************************************************************
#ifndef Exception_h
#define Exception_h

namespace Quinoa {

//! Exception types
enum ExceptionType { CUMULATIVE=0,  //!< Only several will produce a warning
                     WARNING,       //!< Warning: output message
                     ERROR,         //!< Error: output but will not interrupt
                     UNCAUGHT,      //!< Uncaught: will interrupt
                     FATAL          //!< Fatal error: will interrupt
};

//! Error codes for the OS (or whatever calls Quinoa)
enum ErrorCode { NO_ERROR=0,        //!< Everything went fine
                 NONFATAL_ERROR,    //!< Exception occurred but continue
                 FATAL_ERROR,       //!< Fatal error occurred
                 NUM_OF_ERRORS
};

class Driver;

//! Exception base class
class Exception {

  public:
    //! Constructor
    Exception(ExceptionType exception) : m_exception(exception) {}

    //! Copy operator
    Exception(const Exception&);

    //! Destructor
    ~Exception() {}

    //! Handle Exception
    ErrorCode handleException(Driver* driver);

  private:
    //! Don't permit assigment operator
    Exception& operator=(const Exception&);

    //! Exception (CUMULATIVE, WARNING, ERROR, etc.)
    ExceptionType m_exception;
};

} // namespace Quinoa

#endif // Exception_h
