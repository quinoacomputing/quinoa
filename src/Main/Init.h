//******************************************************************************
/*!
  \file      src/Main/Init.h
  \author    J. Bakosi
  \date      Thu 03 Jul 2014 06:23:54 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Common initialization for all mains
  \details   Common initialization for all mains
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

//! Executable types for which an ascii logo is available in Print
enum class HeaderType : uint8_t { QUINOA=0,
                                  RNGTEST,
                                  MESHCONV };

//! Echo program title
void echoHeader( const tk::Print& print, HeaderType header );

//! Echo build environment
void echoBuildEnv( const tk::Print& print,
                   const std::string& executable,
                   void (*echoTPL)(const Print& print) );

//! Echo runtime environment
void echoRunEnv( const tk::Print& print, int argc, char** argv );

//! Generic Main() used for all executables for code-reuse and a consistent look
template< class Driver >
Driver Main( int argc, char* argv[],
             HeaderType header,
             const std::string& executable,
             void (*echoTPL)(const Print&) = [](const Print&){} )
{
  try {

    // Install our own new-handler
    std::set_new_handler( newHandler );
    // Install our own terminate-handler
    std::set_terminate( terminateHandler );
    // Install our own unexpected-handler
    std::set_unexpected( unexpectedHandler );

    // Create pretty printer
    Print print;

    // Echo program header
    echoHeader( print, header );

    // Echo environment
    print.part( "Environment" );
    echoBuildEnv( print, executable, echoTPL ); //!< Build environment
    echoRunEnv( print, argc, argv );            //!< Runtime environment

  } catch (...) { processException(); }

  // Create and return driver
  return Driver( argc, argv );
}

} // tk::

#endif // Init_h
