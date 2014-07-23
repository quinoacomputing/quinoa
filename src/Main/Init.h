//******************************************************************************
/*!
  \file      src/Main/Init.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:23:14 AM MDT
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
#include <tkTags.h>

namespace tk {

//! Wrapper for POSIX API's getcwd() from unistd.h
std::string workdir();

//! Wrapper for the standard C library's gettimeofday() from
std::string curtime();

//! Executable types for which an ascii logo is available in Print
enum class HeaderType : uint8_t { QUINOA=0,
                                  RNGTEST,
                                  UNITTEST,
                                  MESHCONV };

//! Echo program title
void echoHeader( const tk::Print& print, HeaderType header );

//! Echo build environment
void echoBuildEnv( const tk::Print& print,
                   const std::string& executable,
                   void (*echoTPL)(const Print& print) );

//! Echo runtime environment
void echoRunEnv( const tk::Print& print, int argc, char** argv, bool verbose );

//! Generic Main() used for all executables for code-reuse and a uniform output
template< class Driver, class Printer, class CmdLine >
Driver Main( int argc, char* argv[],
             const CmdLine& cmdline,
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
    // Build environment
    echoBuildEnv( print, executable, echoTPL );
    // Runtime environment
    echoRunEnv( print, argc, argv, cmdline.template get< tag::verbose >() );

  } catch (...) { processException(); }

  // Create and return driver
  return Driver( print, cmdline );
}

} // tk::

#endif // Init_h
