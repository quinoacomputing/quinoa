// *****************************************************************************
/*!
  \file      src/Base/Exception.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Exception class declaration
  \details   Exception class declaration. The basic functionality provided by
    the Exception class is to facilitate printing out a message, together with
    the location of the exception (file, line, funcion name), as well as a call
    trace if available, when an exception is thrown. This file also defines
    three macros, Throw, Assert, and ErrChk, that help simplifying client code
    throwing exceptions.
*/
// *****************************************************************************
#ifndef Exception_h
#define Exception_h

#include <exception>
#include <cstdlib>
#include <string>

//! Toolkit declarations and definitions for general purpose utilities
namespace tk {

//! \brief Throw macro that always throws an exception
//! \details Throw Exception with arguments passed in. Add source filename,
//!   function name, and line number where exception occurred. This macro
//!   facilitates a throw of Exception that is somehwat cleaner at the point
//!   of invocation than a direct throw of Exception, as it hides the
//!   file:func:line arguments. Whenever is possible, it should be used via the
//!   Assert and ErrChk macros defined below.
#define Throw(...) \
   throw tk::Exception(__VA_ARGS__, __FILE__, __PRETTY_FUNCTION__, __LINE__)

//! \brief Assert macro that only throws an exception if expr fails.
//! \details If NDEBUG is defined (e.g. cmake's RELEASE or OPTIMIZED mode), do
//!    nothing, expr is not evaluated. If NDEBUG is not defined, evaluate expr.
//!    If expr is true, do nothing. If expr is false, throw Exception with
//!    arguments passed in. The behavior is similar to libc's assert macro, but
//!    throwing an Exception instead will also generate a nice call-trace and
//!    will attempt to free memory. This macro should be used to detect
//!    programmer errors.
#ifdef NDEBUG
#  define Assert(expr, ...) (static_cast<void>(0))
#else  // NDEBUG
#  define Assert(expr, ...) \
   ((expr) ? static_cast<void>(0) : Throw(__VA_ARGS__))
#endif // NDEBUG

//! \brief ErrChk macro that only throws an exception if expr fails.
//! \details The behavior of this macro is the same whether NDEBUG is defined or
//!    not: expr is always evaluated. If expr is true, do nothing. If expr is
//!    false, throw Exception with arguments passed in. This macro should be
//!    used to detect user or runtime errors.
#define ErrChk(expr, ...) \
   ((expr) ? static_cast<void>(0) : Throw(__VA_ARGS__))

//! Error codes for the OS (or whatever calls us)
enum ErrCode { SUCCESS = EXIT_SUCCESS, //!< Everything went fine
               FAILURE = EXIT_FAILURE  //!< Exceptions occurred
};

//! \brief Basic exception class for producing file:func:line info + call trace
//! \details The basic functionality provided by the Exception class is to
//!   facilitate printing out a message, together with the location of the
//!   exception (file, line, funcion name), as well as a call trace if
//!   available, when an exception is thrown.
class Exception : public std::exception {

  public:
    //! Constructor
    explicit Exception( std::string&& message,
                        std::string&& file = "",
                        std::string&& function = "",
                        unsigned int line = 0 ) noexcept;

    //! Destructor
    virtual ~Exception() noexcept;

    //! Force move constructor for throws
    Exception(Exception&&) = default;

    //! Redefine std::exception's what()
    //! \return C-style string to exception message
    virtual const char* what() const noexcept { return m_message.c_str(); }

    //! Handle Exception
    virtual ErrCode handleException() noexcept;

    //! Accessor to function name
    //! \return Reference to function name in which the exception occurred
    const std::string& func() const noexcept { return m_func; }

  private:
    // Use move constructor by default
    //! Don't permit copy constructor
    Exception(const Exception&) = delete;
    //! Don't permit copy assignment
    Exception& operator=(const Exception&) = delete;
    //! Don't permit move assignment
    Exception& operator=(Exception&&) = delete;

    //! Save call trace
    void saveTrace() noexcept;

    //! Demangle and Echo call trace
    void echoTrace() noexcept;

    const std::string m_file;   //!< Source file where exception is occurred
    const std::string m_func;   //!< Function name where exception is occurred
    const unsigned int m_line;  //!< Source line where exception is occurred

    std::string m_message;      //!< Error message
    void* m_addrList[128];      //!< Call-stack before exception
    int m_addrLength;           //!< Number of stack frames
    char** m_symbolList;        //!< Symbol list of stack entries
};

} // tk::

#endif // Exception_h
