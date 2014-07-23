//******************************************************************************
/*!
  \file      src/Main/Init.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:25:07 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Common initialization for mains
  \details   Common initialization for mains
*/
//******************************************************************************

#include <sstream>
#include <ctime>

#include <unistd.h>

#include <Config.h>

#ifdef HAS_BOOST_SYSTEM   // Boost.Asio requires the Boost.System library
#include <boost/asio/ip/host_name.hpp>
#endif

#include <Init.h>
#include <Exception.h>

std::string tk::workdir()
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

std::string tk::curtime()
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

void tk::echoHeader( const Print& print, HeaderType header )
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

void tk::echoBuildEnv( const Print& print,
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

void tk::echoRunEnv( const Print& print, int argc, char** argv, bool verbose )
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
