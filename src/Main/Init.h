//******************************************************************************
/*!
  \file      src/Main/Init.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 03:41:24 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Common initialization for all mains
  \details   Common initialization for all mains.
             Note that including this header requires the upstream definition of
             mainProxy, i.e., it must be included after a Charm++-generated
             decl.h that declares the main charm module unique to an executable.
             This is not very elegant, but this way the code below is reused
             among different executables.
*/
//******************************************************************************
#ifndef Init_h
#define Init_h

#include <string>
#include <sstream>
#include <ctime>
#include <unistd.h>

#ifdef HAS_BOOST_SYSTEM   // Boost.Asio requires the Boost.System library
#include <boost/asio/ip/host_name.hpp>
#endif

#include <Config.h>
#include <Exception.h>
#include <Print.h>
#include <tkTags.h>
#include <Handler.h>

namespace tk {

//! Executable types for which an ascii logo is available in Print
enum class HeaderType : uint8_t { QUINOA=0,
                                  RNGTEST,
                                  UNITTEST,
                                  MESHCONV };


static std::string workdir()
//******************************************************************************
//  Wrapper for POSIX API's getcwd() from unistd.h
//! \author  J. Bakosi
//******************************************************************************
{
  char cwd[1024];

  if (getcwd(cwd, sizeof(cwd)) != NULL)
    return std::string(cwd);
  else
    Throw( "Error from POSIX API's getcwd()" );
}

static std::string curtime()
//******************************************************************************
//  Wrapper for the standard C library's gettimeofday() from
//! \author  J. Bakosi
//******************************************************************************
{
  time_t current_time;
  char* c_time_string;

  // Obtain current time as seconds elapsed since the Epoch
  current_time = time(NULL);

  if (current_time == ((time_t)-1))
    Throw( "Failure to compute the current time" );

  // Convert to local time format
  c_time_string = ctime(&current_time);

  if (c_time_string == NULL)
    Throw( "Failure to convert the current time" );

  // Convert to std::string and remove trailing newline
  std::string str(c_time_string);
  str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());

  return str;
}

static void echoHeader( const Print& print, HeaderType header )
//******************************************************************************
//  Echo program header
//! \author  J. Bakosi
//******************************************************************************
{
  if ( header == HeaderType::QUINOA )
    print.headerQuinoa();
  else if ( header == HeaderType::RNGTEST )
    print.headerRNGTest();
  else if ( header == HeaderType::UNITTEST )
    print.headerUnitTest();
  else if ( header == HeaderType::MESHCONV )
    print.headerMeshConv();
  else
    Throw( "Header not available" );
}

static void echoBuildEnv( const Print& print,
                          const std::string& executable,
                          void (*echoTPL)(const Print& print) )
//******************************************************************************
//  Echo build environment
//! \details Echo information read from [build]/Base/Config.h filled by
//!          CMake based on src/Main/Config.h.in
//! \author  J. Bakosi
//******************************************************************************
{
  print.section( "Build environment" );
  print.item( "Hostname", BUILD_HOSTNAME );
  print.item( "Executable", executable );
  print.item( "Version", VERSION );
  print.item( "Release", RELEASE );
  print.item( "Revision", GIT_COMMIT );
  print.item( "CMake build type", BUILD_TYPE );

#ifdef NDEBUG
  print.item( "Asserts", "off (turn on: CMAKE_BUILD_TYPE=DEBUG)" );
  print.item( "Exception trace", "off (turn on: CMAKE_BUILD_TYPE=DEBUG)" );
#else
  print.item( "Asserts", "on (turn off: CMAKE_BUILD_TYPE=RELEASE)" );
  print.item( "Exception trace", "on (turn off: CMAKE_BUILD_TYPE=RELEASE)" );
#endif

  print.item( "MPI C++ wrapper", MPI_COMPILER );
  print.item( "Underlying C++ compiler", COMPILER );
  print.item( "Build date", BUILD_DATE );

  // Echo info on TPLs used
  echoTPL( print );
}

static void echoRunEnv( const Print& print, int argc, char** argv, bool verbose )
//******************************************************************************
//  Echo runtime environment
//! \author  J. Bakosi
//******************************************************************************
{
  print.section( "Run-time environment" );

  #ifdef HAS_BOOST_SYSTEM
  print.item( "Hostname", boost::asio::ip::host_name() );
  #endif

  print.item( "Date, time", curtime() );
  print.item( "Work directory", workdir() );
  print.item( "Executable (rel. to work dir)", argv[0] );

  print.item("Command line arguments" );
  print << '\'';
  if (argc>1) {
    for (auto i=1; i<argc-1; ++i) {
      print << std::string( argv[i] ) + ' ';
    }
    print << std::string( argv[argc-1] );
  }
  print << "'\n";

  print.item( "Output", verbose ? "verbose (quiet: omit -v)" : "quiet" );
}

static void newHandler()
//******************************************************************************
// Quinoa's own new handler
//! \details This function is automatically called by allocation functions
//! whenever a memory allocation attempt fails. We simply throw an Exception
//! (derived from std::exception, that produces a trace. Exception safety:
//! no-throw guarantee: this function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {
    Throw( "Cannot allocate memory" );
  } catch (...) { processException(); }
}

static void terminateHandler()
//******************************************************************************
// Quinoa's own terminate handler
//! \details This function is automatically called when the exception handling
//! process has to be abandoned for some reason. This happens when no catch
//! handler can be found for a thrown exception, or for some other exceptional
//! circumstance that makes impossible to continue the exception handling
//! process.We simply throw an Exception (derived from std::exception, that
//! produces a trace. Exception safety: no-throw guarantee: this function never
//! throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {
    Throw( "Terminate called by uncaught exception" );
  } catch (...) { processException(); }
}

static void unexpectedHandler()
//******************************************************************************
// Quinoa's own unexpected handler
//! \details This function is automatically called when a function throws an
//! exception that is not in its dynamic-exception-specification (i.e., in its
//! throw specifier). The use of dynamic-exception-specifiers is deprecated
//! (since C++11). We simply throw an Exception (derived from std::exception,
//! that produces a trace. Exception safety: no-throw guarantee: this function
//! never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  try {
    Throw( "Unexpected exception. Exception specifiers are deprecated since "
           "C++11" );
  } catch (...) { processException(); }
}

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
