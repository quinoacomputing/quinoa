//******************************************************************************
/*!
  \file      src/Main/Init.h
  \author    J. Bakosi
  \date      Sun 06 Jul 2014 08:14:54 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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
template< class Driver, class Printer >
Driver Main( int argc, char* argv[],
             HeaderType header,
             const std::string& executable,
             const Printer& print,
             void (*echoTPL)(const Print&) = [](const Print&){} )
{
  try {

    // Install our own new-handler
    std::set_new_handler( newHandler );
    // Install our own terminate-handler
    std::set_terminate( terminateHandler );
    // Install our own unexpected-handler
    std::set_unexpected( unexpectedHandler );

    // Echo program header
    echoHeader( print, header );

    // Echo environment
    print.part( "Environment" );
    echoBuildEnv( print, executable, echoTPL ); //!< Build environment
    echoRunEnv( print, argc, argv );            //!< Runtime environment

  } catch (...) { processException(); }

  // Create and return driver
  return Driver( argc, argv, print );
}

} // tk::

#endif // Init_h
