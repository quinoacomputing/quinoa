//******************************************************************************
/*!
  \file      src/Main/Init.C
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 02:36:30 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Common initialization for mains
  \details   Common initialization for mains
*/
//******************************************************************************

#include <ctime>
#include <unistd.h>

#include <Init.h>
#include <Exception.h>
#include <Config.h>

//! Common initialization routines
namespace init {

using namespace quinoa;

std::string workdir()
//******************************************************************************
//  Wrapper for POSIX API's getcwd() from unistd.h
//! \author  J. Bakosi
//******************************************************************************
{
  char cwd[1024];

  if (getcwd(cwd, sizeof(cwd)) != NULL) {
    return std::string(cwd);
  } else {
    Throw(ExceptType::WARNING, std::string("Error from POSIX API's getcwd()"));
  }
}

std::string curtime()
//******************************************************************************
//  Wrapper for the standard C library's gettimeofday() from
//! \author  J. Bakosi
//******************************************************************************
{
  time_t current_time;
  char* c_time_string;

  // Obtain current time as seconds elapsed since the Epoch
  current_time = time(NULL);

  if (current_time == ((time_t)-1)) {
    Throw(ExceptType::WARNING, "Failure to compute the current time.");
  }

  // Convert to local time format
  c_time_string = ctime(&current_time);

 if (c_time_string == NULL) {
   Throw(ExceptType::WARNING, "Failure to convert the current time.");
 }

 // Convert to std::string and remove trailing newline
 std::string str(c_time_string);
 str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());

 return str;
}

void echoHeader(const Print& print, const std::string& title)
//******************************************************************************
//  Echo program title
//! \author  J. Bakosi
//******************************************************************************
{
  print.header(title);
}

void echoBuildEnv(const Print& print, const std::string& executable)
//******************************************************************************
//  Echo build environment
//! \details Echo information read from [build]/Base/Config.h filled by
//!          CMake based on src/Main/Config.h.in
//! \author  J. Bakosi
//******************************************************************************
{
  print.section("Build environment");
  print.item("Executable", executable);
  print.item("Version", VERSION);
  print.item("Release", RELEASE);
  print.item("Revision", GIT_COMMIT);
  print.item("CMake build type", BUILD_TYPE);
#ifdef NDEBUG
  print.item("Asserts", "off (CMAKE_BUILD_TYPE=DEBUG turns this on)");
  print.item("Exception trace", "off (CMAKE_BUILD_TYPE=DEBUG turns this on)");
#else
  print.item("Asserts", "on (CMAKE_BUILD_TYPE=RELEASE turns this off)");
  print.item("Exception trace", "on (CMAKE_BUILD_TYPE=RELEASE turns this off)");
#endif
  print.item("MPI C++ wrapper", MPI_COMPILER);
  print.item("Underlying C++ compiler", COMPILER);
  print.item("Build date", BUILD_DATE);
}

void echoRunEnv(const Print& print, int argc, char** argv)
//******************************************************************************
//  Echo runtime environment
//! \author  J. Bakosi
//******************************************************************************
{
  print.section("Run-time environment");
  print.item("Date, time", curtime());
  print.item("Work directory", workdir());
  print.item("Executable (rel. to work dir)", argv[0]);

  print.item("Command line arguments");
  print.raw('\'');
  if (argc>1) {
    for (int i=1; i<argc-1; ++i) {
      print.raw(std::string(argv[i]) + ' ');
    }
    print.raw(std::string(argv[argc-1]));
  }
  print.raw("'\n");
}

} // init::
