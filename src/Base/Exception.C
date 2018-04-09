// *****************************************************************************
/*!
  \file      src/Base/Exception.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Exception class definition
  \details   Exception class definition
*/
// *****************************************************************************

#include <string>
#include <sstream>
#include <type_traits>
#include <cstdio>
#include <cstddef>
#include <cxxabi.h>
#include <execinfo.h>

#include "QuinoaConfig.h"
#include "Exception.h"

using tk::Exception;

Exception::Exception( std::string&& message,
                      std::string&& file,
                      std::string&& function,
                      unsigned int line ) noexcept
// *****************************************************************************
//  Constructor: generate error message
//! \param[in] message String (moved from) with an error message
//! \param[in] file String (moved from) with the file name in which the
//!   exception ocurred
//! \param[in] function String (moved from) containing the name of the function
//    in which the exception ocurred
//! \param[in] line Source code line number at which the exception ocurred
//! \details While throwing exceptions is possible calling this constructor, the
//!   easiest and recommend way is to use the Assert, ErrChk, and Throw macros.
//!   Exception safety: no-throw guarantee: this member function never throws
//!   exceptions.
//! \see Assert, ErrChk, Throw
// *****************************************************************************
try :
  m_file( std::move(file) ),
  m_func( std::move(function) ),
  m_line( std::move(line) ),
  m_message( std::move(message) ),
  m_addrLength( 0 ),
  m_symbolList( nullptr )
{

  // Construct exception message
  std::stringstream s;
  s << m_message << std::endl;
  if (line) {
    s << ">>> Exception at " << m_file << ":" << m_line << ": " << m_func;
  } else {
    s << ">>> No file:line:func information from exception";
  }
  m_message = s.str();

  // Save call-trace
  saveTrace();

} // Catch std::exception
  catch (exception& se) {
    // Emit warning and continue
    fprintf( stderr, "RUNTIME ERROR in Exception constructor: %s\n", se.what() );
  }
  // Catch uncaught exceptions
  catch (...) {
    // Emit warning and continue
    fprintf( stderr, "UNKNOWN EXCEPTION in Exception constructor\n" );
  }

Exception::~Exception() noexcept
// *****************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: this member function never
//!   throws exceptions.
// *****************************************************************************
{
  // allocated by execinfo.h's backtrace_symbols() in Exception::saveTrace()
  free(m_symbolList);
}

void
Exception::saveTrace() noexcept
// *****************************************************************************
//  Save call-trace
//! \details Exception safety: no-throw guarantee: this member function never
//!   throws exceptions. For more information see the libc manual at
//!   http://www.gnu.org/software/libc/manual/html_node/Backtraces.html.
//!   Requires stdio.h, execinfo.h.
// *****************************************************************************
{
#ifndef HOST_OS_ALPINE
  // Retrieve current stack addresses
  m_addrLength = backtrace(m_addrList, sizeof(m_addrList)/sizeof(void*));
#endif

  // Resolve addresses into strings containing "filename(function+address)"
  if (m_addrLength > 0)
    m_symbolList = backtrace_symbols(m_addrList, m_addrLength);
}

void
Exception::echoTrace() noexcept
// *****************************************************************************
//  Demangle and echo call trace
//! \details Exception safety: no-throw guarantee: this member function never
//!   throws exceptions. Credit goes to Timo Bingmann at http://idlebox.net,
//!   published under the WTFPL v2.0. For more information see
//!   * http://panthema.net/2008/0901-stacktrace-demangled
//!   * http://panthema.net/2008/0901-stacktrace-demangled/cxa_demangle.html
//!   * http://gcc.gnu.org/onlinedocs/libstdc++/latest-doxygen
// *****************************************************************************
{
  // Allocate string which will be filled with the demangled function name
  size_t funcnamesize = 256;
  char* funcname = static_cast< char* >( malloc(funcnamesize) );

  // Iterate over the returned symbol lines. skip the first two, these are the
  // addresses of Exception::saveTrace() and the Exception constructor
  for (int i=2; i<m_addrLength; ++i) {
    char *begin_name = nullptr, *begin_offset = nullptr, *end_offset = nullptr;

    // Find parentheses and +address offset surrounding the mangled name:
    // ./module(function+0x15c) [0x8048a6d]
    for (char *p = m_symbolList[i]; *p; ++p) {
      if (*p == '(')
        begin_name = p;
      else if (*p == '+')
        begin_offset = p;
      else if (*p == ')' && begin_offset) {
        end_offset = p;
        break;
      }
    }

    if (begin_name && begin_offset && end_offset && begin_name < begin_offset) {
      *begin_name++ = '\0';
      *begin_offset++ = '\0';
      *end_offset = '\0';

      // Mangled name is now in [begin_name, begin_offset) and caller
      // offset in [begin_offset, end_offset). now apply __cxa_demangle()
      int status;
      char* ret = abi::__cxa_demangle(begin_name,
                                      funcname, &funcnamesize, &status);

      if (status == 0) {
        funcname = ret; // use possibly realloc()-ed string
        fprintf( stderr, ">>> %s : %s+%s\n", m_symbolList[i], funcname,
                         begin_offset );
      } else {
        // Demangling failed. Output function name as a C function with no
        // arguments
        fprintf( stderr, ">>> %s : %s()+%s\n", m_symbolList[i], begin_name,
                         begin_offset);
      }
    } else {
      // Couldn't parse the line? Print the whole line
      fprintf( stderr, ">>> %s\n", m_symbolList[i] );
    }
  }

  free(funcname);
}

tk::ErrCode
Exception::handleException() noexcept
// *****************************************************************************
//  Handle Exception: Print cumulative message
//! \return Error code, as defined in stdlib.h, i.e., cstdlib
//! \details Exception safety: no-throw guarantee: this member function never
//!   throws exceptions.
// *****************************************************************************
{
  if (m_addrLength > 0) {
    fprintf( stderr, ">>>\n>>> =========== CALL TRACE ===========\n>>>\n" );
    echoTrace();
    fprintf( stderr, ">>>\n>>> ======= END OF CALL TRACE ========\n>>>\n" );
  }
 
  return tk::ErrCode::FAILURE;
}
