//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Thu Sep 19 17:29:02 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>

#include <QuinoaConfig.h>
#include <Base.h>
#include <Handler.h>
#include <Paradigm.h>
#include <Memory.h>
#include <QuinoaDriver.h>
#include <QuinoaPrint.h>

using namespace quinoa;

//! Everything that contributes to the quina executable
namespace quinoa {

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

static void echoRunEnv(const QuinoaPrint& print)
//******************************************************************************
//  Echo runtime environment
//! \author  J. Bakosi
//******************************************************************************
{
  print.section("Run-time environment");
  print.item("Date, time", "...");
  print.item("Working directory", "...");
  print.item("Executable full path", "...");
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
  try {

    // Install our own new-handler
    std::set_new_handler(newHandler);
    // Install our own terminate-handler
    std::set_terminate(terminateHandler);
    // Install our own unexpected-handler
    std::set_unexpected(unexpectedHandler);

    // Create the essentials
    QuinoaPrint print;                  //!< Pretty printer
    Paradigm paradigm(print);           //!< Parallel compute environment
    Memory memory(&paradigm);           //!< Memory manager
    QuinoaControl control;              //!< Controller
    Timer timer;                        //!< Timer

    // Bundle up essentials
    Base base(print, paradigm, memory, control, timer);

    // Echo program name
    echoHeader(print);

    // Echo environment
    print.part("Environment");
    echoBuildEnv(print);                //!< Build environment
    echoRunEnv(print);                  //!< Runtime environment
    paradigm.echo();                    //!< Parallel compute enviroment
    print.endpart();

    // Create driver
    QuinoaDriver driver(argc, argv, base);

    // Execute
    driver.execute();

  } catch (...) {
      processException();
    }

  // Return error code success
  return static_cast<int>(ErrCode::SUCCESS);
}
