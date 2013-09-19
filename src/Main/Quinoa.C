//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Thu Sep 19 10:38:35 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>

#include <QuinoaConfig.h>
#include <Base.h>
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

} // namespace quinoa

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa main
//! \details Quinoa main -- This is where everything starts...
//! \author  J. Bakosi
//******************************************************************************
{
  ErrCode error = ErrCode::SUCCESS;
  try {

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

    // Solve
    driver.execute();

  } // Catch and handle Quina::Exceptions
    catch (Exception& qe) {
      error = qe.handleException();
    }
    // Catch std::exceptions and transform them into Quinoa::Exceptions without
    // file:line:func information
    catch (std::exception& se) {
      Exception qe(ExceptType::RUNTIME, se.what());
      error = qe.handleException();
    }
    // Catch uncaught exceptions and still do cleanup
    catch (...) {
      Exception qe(ExceptType::UNCAUGHT, "Non-standard exception");
      error = qe.handleException();
    }

  // Return error code
  return static_cast<int>(error);
}
