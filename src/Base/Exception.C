//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Wed May 29 09:06:06 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class definition
  \details   Exception base class definition
*/
//******************************************************************************

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cxxabi.h>

#include <execinfo.h>

#include <Driver.h>
#include <Exception.h>

using namespace Quinoa;

Exception::Exception(const ExceptType except,
                     const string& message,
                     const string& file,
                     const string& func,
                     const unsigned int line) noexcept
//******************************************************************************
//  Constructor: generate error message
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
try :
  m_except(move(except)),
  m_file(move(file)),
  m_func(move(func)),
  m_line(move(line)),
  m_message(move(message)),
  m_addrLength(0),
  m_symbolList(nullptr)
{

  // Construct exception message
  stringstream s;
  s << m_message << endl;
  if (line) {
    s << ">>> Exception in " << m_file << ":" << m_line << ": " << m_func;
  } else {
    s << ">>> No file:line:func information from exception";
  }
  m_message = s.str();

  // Save call-trace
  saveTrace();

} // Catch std::exception
  catch (exception& se) {
    // Emit warning and continue
    cout << "RUNTIME ERROR in Exception constructor: " << se.what() << endl;
  }
  // Catch uncaught exceptions
  catch (...) {
    // Emit warning and continue
    cout << "UNKNOWN EXCEPTION in Exception constructor" << endl;
  }

Exception::~Exception() noexcept
//******************************************************************************
//  Destructor
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  // allocated by execinfo.h's backtrace_symbols() in Exception::saveTrace()
  free(m_symbolList);
}

void
Exception::saveTrace() noexcept
//******************************************************************************
//  Save call-trace
//! \details No-throw guarantee: this member function never throws exceptions.
//!          For more information see the libc manual at
//!          http://www.gnu.org/software/libc/manual/html_node/Backtraces.html.
//!          Requires stdio.h, execinfo.h.
//! \author J. Bakosi
//******************************************************************************
{
  // Retrieve current stack addresses
  m_addrLength = backtrace(m_addrList, sizeof(m_addrList)/sizeof(void*));

  if (m_addrLength == 0) {
     printf(">>> Call stack is empty, possibly corrupt.\n");
  }

  // Resolve addresses into strings containing "filename(function+address)"
  m_symbolList = backtrace_symbols(m_addrList, m_addrLength);
}

void
Exception::echoSymbols() noexcept
//******************************************************************************
//  Echo call trace as symbol list
//! \details No-throw guarantee: this member function never throws exceptions.
//!          For more information see the libc manual at
//!          http://www.gnu.org/software/libc/manual/html_node/Backtraces.html
//! \author J. Bakosi
//******************************************************************************
{
  // Echo trace
  for (int i=0; i<m_addrLength; ++i) {
    printf("%s\n", m_symbolList[i]);
  }
}

void
Exception::echoTrace() noexcept
//******************************************************************************
//  Demangle and echo call trace
//! \details No-throw guarantee: this member function never throws exceptions.
//!   Credit goes to Timo Bingmann at http://idlebox.net, published under the
//!   WTFPL v2.0. For more information see
//!   http://panthema.net/2008/0901-stacktrace-demangled 
//!   http://panthema.net/2008/0901-stacktrace-demangled/cxa_demangle.html
//!   http://gcc.gnu.org/onlinedocs/libstdc++/latest-doxygen
//! \author J. Bakosi
//******************************************************************************
{
  // Allocate string which will be filled with the demangled function name
  size_t funcnamesize = 256;
  char* funcname = (char*)malloc(funcnamesize);

  // Iterate over the returned symbol lines. skip the first two, these are the
  // addresses of Exception::saveTrace() and the Exception constructor
  for (int i=2; i<m_addrLength; ++i) {
    char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

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
        printf(">>> %s : %s+%s\n", m_symbolList[i], funcname, begin_offset);
      } else {
        // Demangling failed. Output function name as a C function with no
        // arguments
        printf(">>> %s : %s()+%s\n", m_symbolList[i], begin_name, begin_offset);
      }
    } else {
      // Couldn't parse the line? Print the whole line
      printf(">>> %s\n", m_symbolList[i]);
    }
  }

  free(funcname);
}

void
Exception::echo(const char* msg) noexcept
//******************************************************************************
//  Echo message
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  printf(">>> %s: %s\n>>> CALL TRACE: ======================================="
         "==============\n", msg, what());
  echoTrace();
  printf(">>> ==============================================================="
         "==\n");
}

ErrCode
Exception::handleException(Driver* const driver) noexcept
//******************************************************************************
//  Handle Exception: Print cumulative message and handle criticality
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  if (m_except == ExceptType::WARNING) {

    echo("WARNING");
    return ErrCode::NONFATAL;

  } else if (m_except == ExceptType::CUMULATIVE) {

    echo("CUMULATIVE ERROR");
    return ErrCode::NONFATAL;

  } else if (m_except == ExceptType::ERROR) {

    echo("ERROR");
    return ErrCode::NONFATAL;

  } else if (m_except == ExceptType::FATAL) {

    echo("FATAL ERROR");
    printf(">>> Attempting cleanup & graceful exit...\n");
    if (driver) driver->finalize();
    return ErrCode::FATAL_ERROR;

  } else if (m_except == ExceptType::RUNTIME) {

    echo("RUNTIME ERROR");
    printf(">>> Attempting cleanup & graceful exit...\n");
    if (driver) driver->finalize();
    return ErrCode::FATAL_ERROR;

  } else if (m_except == ExceptType::UNCAUGHT) {

    echo("UNKNOWN ERROR");
    printf(">>> Attempting cleanup & graceful exit...\n");
    if (driver) driver->finalize();
    return ErrCode::FATAL_ERROR;

  } else {
    return ErrCode::NONFATAL;
  }
}
