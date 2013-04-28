//******************************************************************************
/*!
  \file      src/Base/Exception.h
  \author    J. Bakosi
  \date      Sat 27 Apr 2013 08:33:05 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class declaration
  \details   Exception base class declaration
*/
//******************************************************************************
#ifndef Exception_h
#define Exception_h

#include <QuinoaTypes.h>

namespace Quinoa {

// Throw macro that always throws an exception:
// If NDEBUG is defined (e.g. RELEASE/OPTIMIZED mode), do nothing.  If NDEBUG
// is not defined, throw exception passed in as argument 2 with exception-
// arguments passed in as arguments 2+. Add source filename, function name, and
// line number where exception occurred.
#ifdef NDEBUG
#  define Throw(exception, ...) (static_cast<void>(0))
#else  // NDEBUG
#  define Throw(exception, ...) \
   throw exception(__VA_ARGS__, __FILE__, __PRETTY_FUNCTION__, __LINE__)
#endif // NDEBUG

// Assert macro that only throws an exception if expr fails:
// If NDEBUG is defined (e.g. RELEASE/OPTIMIZED mode), do nothing, expr is not
// evaluated.  If NDEBUG is not defined, evaluate expr. If expr is true, do
// nothing. If expr is false, throw exception passed in as argument 2 with
// exception-arguments passed in as arguments 2+.
#ifdef NDEBUG
#  define Assert(expr, exception, ...) (static_cast<void>(0))
#else  // NDEBUG
#  define Assert(expr, exception, ...) \
   ((expr) ? static_cast<void>(0) : Throw(exception,__VA_ARGS__))
#endif // NDEBUG

// Errchk macro that only throws an exception if expr fails:
// If NDEBUG is defined (e.g. RELEASE/OPTIMIZED mode), expr is still evaluated.
// If NDEBUG is not defined, evaluate expr. If expr is true, do nothing. If
// expr is false, throw exception passed in as argument 2 with exception-
// arguments passed in as arguments 2+.
#ifdef NDEBUG
#  define Errchk(expr, exception, ...) (expr)
#else  // NDEBUG
#  define Errchk(expr, exception, ...) \
   ((expr) ? static_cast<void>(0) : Throw(exception,__VA_ARGS__))
#endif // NDEBUG

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
class Exception {

  public:
    //! Constructor: empty with defaults
    explicit Exception(const ExceptType except, const string& msg = "") :
      m_file(""),
      m_func(""),
      m_line(0),
      m_message(msg),
      m_except(except) {}

    //! Destructor
    virtual ~Exception() noexcept = default;

    //! Handle Exception passing pointer to driver
    virtual ErrCode handleException(Driver* const driver);

    //! Generate generic exception message
    string genericWhat() const;

    //! Generate std::runtime_error exception message
    string runtimeWhat() const;

  protected:
    //! Constructor: initialized from children
    explicit Exception(const ExceptType except,
                       const string& file,
                       const string& func,
                       const unsigned int& line,
                       const string& message) noexcept :
      m_file(file),
      m_func(func),
      m_line(line),
      m_message(message),
      m_except(except) {}

    //! Permit copy constructor only for children
    Exception(const Exception&) = default;
    //! Permit move constructor only for children
    Exception(Exception&&) = default;

    //! Augment message after its construction
    //! \param[in]  message  Message to add to end of message
    void augment(const string& message) { m_message += message; }

    const string m_file;  //!< Source file where the exception is occurred
    const string m_func;  //!< Functionn name in which the exception is occurred
    const int m_line;     //!< Source line where the exception is occurred

  private:
    //! Don't permit copy assignment
    Exception& operator=(const Exception&) = delete;
    //! Don't permit move assignment
    Exception& operator=(Exception&&) = delete;

    string m_message;     //!< Error message (constructed along the tree)

    //! Exception type (WARNING, CUMULATIVE, ERROR, etc.)
    const ExceptType m_except;
};

} // namespace Quinoa

#endif // Exception_h
