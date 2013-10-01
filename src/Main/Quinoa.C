//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 10:25:01 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>

#include <unistd.h>
#include <time.h>

#include <QuinoaConfig.h>
#include <Base.h>
#include <Handler.h>
#include <Paradigm.h>
#include <Memory.h>
#include <QuinoaDriver.h>
#include <QuinoaPrint.h>
#include <Quinoa/CmdLine/Keywords.h>

using namespace quinoa;

//! Everything that contributes to the quina executable
namespace quinoa {

static std::string workdir()
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

  if (current_time == ((time_t)-1)) {
    Throw(ExceptType::WARNING, "Failure to compute the current time.");
  }

  // Convert to local time format
  c_time_string = ctime(&current_time);

 if (c_time_string == NULL) {
   Throw(ExceptType::FATAL, "Failure to convert the current time.");
 }

 // Convert to std::string and remove trailing newline
 std::string str(c_time_string);
 str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());

 return str;
}

static void echoHeader(const QuinoaPrint& print)
//******************************************************************************
//  Echo Name
//! \author  J. Bakosi
//******************************************************************************
{
  print.header("Quinoa: Lagrangian particle hydrodynamics");
}

static void echoBuildEnv(const QuinoaPrint& print)
//******************************************************************************
//  Echo build environment
//! \details Echo information read from [build]/Base/QuinoaConfig.h filled by
//!          CMake based on src/MainQuinoaConfig.h.in
//! \author  J. Bakosi
//******************************************************************************
{
  print.section("Build environment");
  print.item("Executable", QUINOA_EXECUTABLE);
  print.item("Version", QUINOA_VERSION);
  print.item("Release", QUINOA_RELEASE);
  print.item("Revision", QUINOA_GIT_COMMIT);
  print.item("CMake build type", QUINOA_BUILD_TYPE);
#ifdef NDEBUG
  print.item("Asserts", "off");
#else  // NDEBUG
  print.item("Asserts", "on");
#endif // NDEBUG
  print.item("MPI C++ wrapper", QUINOA_MPI_COMPILER);
  print.item("Underlying C++ compiler", QUINOA_COMPILER);
  print.item("Build date", QUINOA_BUILD_DATE);
}

static void echoRunEnv(const QuinoaPrint& print, const char* executable)
//******************************************************************************
//  Echo runtime environment
//! \author  J. Bakosi
//******************************************************************************
{
  print.section("Run-time environment");
  print.item("Date, time", curtime());
  print.item("Work directory", workdir());
  print.item("Executable (rel. to work dir)", executable);
  print.item("Command line arguments", "...");
}

} // quinoa::

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa main
//! \details Quinoa main -- This is where everything starts...
//! \author  J. Bakosi
//******************************************************************************
{
  ErrCode error = ErrCode::SUCCESS;

  try {

    // Install our own new-handler
    std::set_new_handler(newHandler);
    // Install our own terminate-handler
    std::set_terminate(terminateHandler);
    // Install our own unexpected-handler
    std::set_unexpected(unexpectedHandler);

    // Create the essentials
    InputDeck control;                      //!< Parsed input deck
    QuinoaPrint print(control);             //!< Pretty printer
    Paradigm paradigm(print);               //!< Parallel compute environment
    Memory memory(&paradigm);               //!< Memory manager
    Timer timer;                            //!< Timer

    // Bundle up essentials
    Base base(print, paradigm, memory, control, timer);

    // Echo program name
    echoHeader(print);

    // Echo environment
    print.part("Environment");
    echoBuildEnv(print);                //!< Build environment
    paradigm.echo();                    //!< Parallel compute enviroment
    echoRunEnv(print, argv[0]);         //!< Runtime environment

    // Create driver
    QuinoaDriver driver(argc, argv, base);

    // Execute
    driver.execute();

  } catch (...) {
      error = processException();
    }

  // Return error code success
  return static_cast<int>(error);
}
