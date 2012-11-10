//******************************************************************************
/*!
  \file      src/Base/Exception.h
  \author    J. Bakosi
  \date      Fri 09 Nov 2012 07:14:50 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class declaration
  \details   Exception base class declaration
*/
//******************************************************************************
#ifndef Exception_h
#define Exception_h

#include <iostream>

#include <QuinoaTypes.h>

namespace Quinoa {

//! Exception types
// ICC: no strongly typed enums yet
enum ExceptType { CUMULATIVE=0,  //!< Only several will produce a warning
                  WARNING,       //!< Warning: output message
                  ERROR,         //!< Error: output but will not interrupt
                  FATAL,         //!< Fatal error: will interrupt
                  UNCAUGHT,      //!< Uncaught: will interrupt
                  NUM_EXCEPT
};

//! Error codes for the OS (or whatever calls Quinoa)
enum ErrCode { NO_ERROR=0,       //!< Everything went fine
               NONFATAL,         //!< Exception occurred but continue
               FATAL_ERROR,      //!< Fatal error occurred
               NUM_ERR_CODE
};

class Driver;

//! Exception base class
class Exception {

  public:
    //! Constructor
    Exception(ExceptType except) : m_except(except) {}
    Exception(ExceptType except, const string& msg) :
      m_message(msg), m_except(except) {}

    //! Move constructor, necessary for throws, default compiler generated
    Exception(Exception&&) = default;

    //! Destructor
    virtual ~Exception() {}

    //! Handle Exception passing pointer to driver
    virtual ErrCode handleException(Driver* driver);

  protected:
    //! Don't permit copy constructor
    // ICC: should be deleted and private
    Exception(const Exception&);

    //! Error message (constructed along the inheritance tree)
    string m_message;

  private:
    //! Don't permit copy assignment
    Exception& operator=(const Exception&) = delete;
    //! Don't permit move assignment
    Exception& operator=(Exception&&) = delete;

    //! Exception type (CUMULATIVE, WARNING, ERROR, etc.)
    ExceptType m_except;
};

} // namespace Quinoa

#endif // Exception_h
