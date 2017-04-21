// *****************************************************************************
/*!
  \file      src/Main/Init.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Common initialization routines for main() functions for multiple
     exectuables
  \details   Common initialization routines for main() functions for multiple
     exectuables. The functions in this file are used by multiple execitables
     to ensure code-reuse and a uniform screen-output.
*/
// *****************************************************************************
#ifndef Init_h
#define Init_h

#include <string>
#include <sstream>
#include <ctime>
#include <unistd.h>

#include "QuinoaConfig.h"
#include "Exception.h"
#include "Print.h"
#include "Tags.h"

namespace tk {

//! Executable types for which an ascii logo is available in tk::Print
enum class HeaderType : uint8_t { INCITER=0,
                                  RNGTEST,
                                  UNITTEST,
                                  MESHCONV,
                                  WALKER };


static std::string workdir()
// *****************************************************************************
//! \brief Wrapper for POSIX API's getcwd() from unistd.h
//! \return A stirng containing the current working directory
//! \author J. Bakosi
// *****************************************************************************
{
  char cwd[1024];

  if ( getcwd(cwd, sizeof(cwd)) != NULL )
    return std::string( cwd );
  else
    Throw( "Error from POSIX API's getcwd()" );
}

static std::string curtime()
// *****************************************************************************
//! \brief Wrapper for the standard C library's gettimeofday() from
//! \return A stirng containing the current date and time
//! \author J. Bakosi
// *****************************************************************************
{
  time_t current_time;
  char* c_time_string;

  // Obtain current time as seconds elapsed since the Epoch
  current_time = time( NULL );

  if (current_time == static_cast<time_t>(-1))
    Throw( "Failure to compute the current time" );

  // Convert to local time format
  c_time_string = ctime(&current_time);

  if (c_time_string == NULL)
    Throw( "Failure to convert the current time" );

  // Convert to std::string and remove trailing newline
  std::string str( c_time_string );
  str.erase( std::remove(str.begin(), str.end(), '\n'), str.end() );

  return str;
}

static void echoHeader( const Print& print, HeaderType header )
// *****************************************************************************
//! \brief Echo program header
//! \param[in] print Pretty printer
//! \param[in] header Header type enum indicating which header to print
//! \author  J. Bakosi
// *****************************************************************************
{
  if ( header == HeaderType::INCITER )
    print.headerInciter();
  else if ( header == HeaderType::RNGTEST )
    print.headerRNGTest();
  else if ( header == HeaderType::UNITTEST )
    print.headerUnitTest();
  else if ( header == HeaderType::MESHCONV )
    print.headerMeshConv();
  else if ( header == HeaderType::WALKER )
    print.headerWalker();
  else
    Throw( "Header not available" );
}

static void echoBuildEnv( const Print& print, const std::string& executable )
// *****************************************************************************
//! \brief Echo build environment
//! \details Echo information read from <build>/Base/Config.h filled by
//!    CMake based on src/Main/Config.h.in.
//! \param[in] print Pretty printer
//! \param[in] executable Name of the executable
//! \author J. Bakosi
// *****************************************************************************
{
  print.section( "Build environment" );
  print.item( "Hostname", BUILD_HOSTNAME );
  print.item( "Executable", executable );
  print.item( "Version", VERSION );
  if (std::string(GIT_COMMIT).find("NOTFOUND") == std::string::npos)
    print.item( "Revision SHA1", GIT_COMMIT );
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
}

static void echoRunEnv( const Print& print, int argc, char** argv, bool verbose )
// *****************************************************************************
//! \brief Echo runtime environment
//! \param[in] print Pretty printer
//! \param[in] argc Number of command-line arguments to executable
//! \param[in] argv C-style string array to command-line arguments to executable
//! \param[in] verbose True for verbose screen-output
//! \author J. Bakosi
// *****************************************************************************
{
  print.section( "Run-time environment" );

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

//! \brief Generic Main() used for all executables for code-reuse and a uniform
//!    output
//! \details The template arguments configure this Main class that is
//!   practically used instead of the usual main(). This allows code-reuse and a
//!   unfirom screen-output. The template arguments are:
//!   - Driver, specializaing the driver type to be created, see tk::Driver
//!   - Printer, specializaing the pretty printer type to use, see tk::Print
//!   - CmdLine, specializing the command line object storing data parsed from
//!     the command line
//! \param[in] argc Number of command-line arguments to executable
//! \param[in] argv C-style string array to command-line arguments to executable
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \param[in] header Header type enum indicating which executable header to
//!   print
//! \param[in] executable Name of the executable
//! \param[in] print Pretty printer to use
//! \return Instantiated driver object which can then be used to execute()
//!   whatever it is intended to drive
//! \author J. Bakosi
template< class Driver, class Printer, class CmdLine >
Driver Main( int argc, char* argv[],
             const CmdLine& cmdline,
             HeaderType header,
             const std::string& executable,
             const Printer& print )
{
  // Echo program header
  echoHeader( print, header );

  // Echo environment
  print.part( "Environment" );
  // Build environment
  echoBuildEnv( print, executable );
  // Runtime environment
  echoRunEnv( print, argc, argv, cmdline.template get< tag::verbose >() );

  // Create and return driver
  return Driver( print, cmdline );
}

} // tk::

#endif // Init_h
