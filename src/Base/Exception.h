//******************************************************************************
/*!
  \file      src/Base/Exception.h
  \author    J. Bakosi
  \date      Tue May  7 12:26:38 2013
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

//! Throw macro that always throws an exception:
//! Throw Exception with arguments passed in. Add source filename, function
//! name, and line number where exception occurred. This macro facilitates a
//! throw of Quinoa::Exception that is somehwat cleaner at the point of
//! invocation than a direct throw of Exception, as it hides the file:func:line
//! arguments. Whenever is possible, it should be used via the Assert and ErrChk
//! macros defined below.
#define Throw(...) \
   throw Exception(__VA_ARGS__, __FILE__, __PRETTY_FUNCTION__, __LINE__)

//! Assert macro that only throws an exception if expr fails:
//! If NDEBUG is defined (e.g. RELEASE/OPTIMIZED mode), do nothing, expr is not
//! evaluated. If NDEBUG is not defined, evaluate expr. If expr is true, do
//! nothing. If expr is false, throw Exception with arguments passed in.
//! The behavior is similar to libc's assert macro, but throwing an Exception
//! instead will also generate a nice call-trace and will attempt to free
//! memory. This macro should be used to detect programmer errors.
#ifdef NDEBUG
#  define Assert(expr, ...) (static_cast<void>(0))
#else  // NDEBUG
#  define Assert(expr, ...) \
   ((expr) ? static_cast<void>(0) : Throw(__VA_ARGS__))
#endif // NDEBUG

//! ErrChk macro that only throws an exception if expr fails:
//! The behavior is the same whether NDEBUG is defined or not: expr is always
//! evaluated. If expr is true, do nothing. If expr is false, throw Exception
//! with arguments passed in. This macro should be used to detect user/runtime
//! errors.
#define ErrChk(expr, ...) \
   ((expr) ? static_cast<void>(0) : Throw(__VA_ARGS__))

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
                       const std::string& message,
                       const std::string& file = "",
                       const std::string& func = "",
                       const unsigned int line = 0) noexcept;

    //! Destructor
    virtual ~Exception() noexcept;

    //! Force move constructor for throws
    Exception(Exception&&) = default;

    //! Redefine std::exception's what()
    virtual const char* what() const noexcept { return m_message.c_str(); }

    //! Handle Exception passing pointer to driver whose finalize() is called
    virtual ErrCode handleException(Driver* const driver = nullptr) noexcept;

    //! Echo message
    void echo(const char* msg) noexcept;

    //! Accessor to function name
    const std::string& func() const noexcept { return m_func; }

  protected:
    //! Force copy constructor for children
    // ICC: in C++11 this should be deleted and private
    Exception(const Exception&) = default;

  private:
    //! Don't permit copy assignment
    Exception& operator=(const Exception&) = delete;
    //! Don't permit move assignment
    Exception& operator=(Exception&&) = delete;

    //! Save call trace
    void saveTrace() noexcept;

    //! Echo call trace as symbols
    void echoSymbols() noexcept;

    //! Demangle and Echo call trace
    void echoTrace() noexcept;

    const ExceptType m_except;  //!< Exception type (WARNING, etc.)
    const std::string m_file;   //!< Source file where exception is occurred
    const std::string m_func;   //!< Function name where exception is occurred
    const unsigned int m_line;  //!< Source line where exception is occurred

    std::string m_message;      //!< Error message
    void* m_addrList[128];      //!< Call-stack before exception
    int m_addrLength;           //!< Number of stack frames
    char** m_symbolList;        //!< Symbol list of stack entries
};

} // namespace Quinoa

#endif // Exception_h
