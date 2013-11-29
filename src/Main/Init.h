//******************************************************************************
/*!
  \file      src/Main/Init.h
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 03:15:21 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Common initialization for mains
  \details   Common initialization for mains
*/
//******************************************************************************
#ifndef Init_h
#define Init_h

#include <string>

#include <Handler.h>
#include <Print.h>

namespace tk {

//! Wrapper for POSIX API's getcwd() from unistd.h
std::string workdir();

//! Wrapper for the standard C library's gettimeofday() from
std::string curtime();

#ifdef HAS_MKL
//! Echo MKL (Intel Math Kernel Library) version information
void echoMKL(const tk::Print& print, const std::string& title);
#endif

//! Echo Boost library version information
void echoBoost(const tk::Print& print, const std::string& title);

//! Echo program title
void echoHeader(const tk::Print& print, const std::string& title);

//! Echo build environment
void echoBuildEnv(const tk::Print& print, const std::string& executable);

//! Echo runtime environment
void echoRunEnv(const tk::Print& print, int argc, char** argv);

//! Main() templated on driver
template< class Driver >
int Main(int argc, char* argv[],
         const std::string& name,
         const std::string& executable)
{
  ErrCode error = ErrCode::SUCCESS;

  try {

    // Install our own new-handler
    std::set_new_handler(newHandler);
    // Install our own terminate-handler
    std::set_terminate(terminateHandler);
    // Install our own unexpected-handler
    std::set_unexpected(unexpectedHandler);

    // Create pretty printer
    Print print;

    // Echo program name
    echoHeader(print, name);

    // Echo environment
    print.part("Environment");
    echoBuildEnv(print, executable);  //!< Build environment
    echoRunEnv(print, argc, argv);    //!< Runtime environment

    // Create driver
    Driver driver(argc, argv, print);

    // Execute
    driver.execute();

  } catch (...) {
      error = processException();
    }

  // Return error code
  return static_cast<int>(error);
}

} // tk::

#endif // Init_h
