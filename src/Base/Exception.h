//******************************************************************************
/*!
  \file      src/Base/Exception.h
  \author    J. Bakosi
  \date      Sat 04 May 2013 07:50:36 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class declaration
  \details   Exception base class declaration
*/
//******************************************************************************
#ifndef Exception_h
#define Exception_h

#include <exception>
#include <string>

#include <QuinoaTypes.h>

namespace Quinoa {

//! Exception types
// ICC: no strongly typed enums yet
enum ExceptType { WARNING=0,     //!< produces a warning
                  CUMULATIVE,    //!< only several will produce a warning
                  ERROR,         //!< produce error but will not interrupt
                  FATAL,         //!< fatal error, will interrupt
                  RUNTIME,       //!< std::runtime_error, will interrupt
                  UNCAUGHT,      //!< uncaught: will interrupt
                  NUM_EXCEPT
};

//! Error codes for the OS (or whatever calls Quinoa)
enum ErrCode { HAPPY=0,          //!< Everything went fine
               NONFATAL,         //!< Non-fatal exceptions occurred
               FATAL_ERROR,      //!< Fatal error occurred, had to terminate
               NUM_ERR_CODE
};

class Driver;

//! Quinoa::Exception base
class Exception : public std::exception {

  public:
    //! Constructor
    explicit Exception(const ExceptType except,
                       const std::string& message = "",
                       int number = 0) noexcept;

    //! Destructor
    virtual ~Exception() noexcept = default;

    //! Force move constructor for throws
    Exception(Exception&&) = default;

    //! Redefine std::exception's what()
    virtual const char* what() const noexcept { return m_message.c_str(); }

    //! Handle Exception passing pointer to driver
    virtual ErrCode handleException(Driver* const driver) noexcept;

    //! Echo message
    void echo(const char* msg) noexcept;

  protected:
    //! Force copy constructor for children
    // ICC: in C++11 this should be deleted and private
    Exception(const Exception&) = default;

  private:
    //! Don't permit copy assignment
    Exception& operator=(const Exception&) = delete;
    //! Don't permit move assignment
    Exception& operator=(Exception&&) = delete;

    //! Echo call trace
    void trace() noexcept;

    const ExceptType m_except; //! Exception type (WARNING, CUMULATIVE, etc.)
    std::string m_message;     //!< Error message
};

} // namespace Quinoa

#endif // Exception_h
