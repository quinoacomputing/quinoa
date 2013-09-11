//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Wed Sep 11 17:06:14 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>

#include <QuinoaConfig.h>
#include <Paradigm.h>
#include <Memory.h>
#include <QuinoaDriver.h>
#include <QuinoaPrinter.h>

using namespace quinoa;

namespace quinoa {

static void echoTitle(const QuinoaPrinter& print)
//******************************************************************************
//  Echo Name
//! \author  J. Bakosi
//******************************************************************************
{
  print.title("Quinoa: Lagrangian particle hydrodynamics");
}

static void echoBuildEnv(const QuinoaPrinter& print)
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
  print.item("MPI C++ compiler", QUINOA_MPI_COMPILER);
  print.item("Underlying C++ compiler", QUINOA_COMPILER);
  print.item("Build date", QUINOA_BUILD_DATE);
#ifdef NDEBUG
  print.item("Asserts", "off");
#else  // NDEBUG
  print.item("Asserts", "on");
#endif // NDEBUG
}

static void echoRunEnv(const QuinoaPrinter& print)
//******************************************************************************
//  Echo runtime environment
//! \author  J. Bakosi
//******************************************************************************
{
  print.section("Run-time environment");
  print.item("Date & Time", "...");
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
  Memory* memory = nullptr;
  QuinoaDriver* driver = nullptr;

  ErrCode error = ErrCode::SUCCESS;
  try {

    // Create pretty printer
    QuinoaPrinter printer;
    // Echo program name
    echoTitle(printer);
    // Echo build environment
    echoBuildEnv(printer);
    // Echo runtime environment
    echoRunEnv(printer);

    // Query, setup, and echo parallel enviroment
    Paradigm paradigm(printer);
    paradigm.echo();

    // Initialize memory manager
    memory = new (std::nothrow) Memory(&paradigm);
    ErrChk(memory != nullptr, ExceptType::FATAL,
           "No memory for a memory manager?");

    // Allocate and initialize driver
    driver = new (std::nothrow)
             QuinoaDriver(argc, argv, memory, &paradigm, printer);
    ErrChk(driver != nullptr, ExceptType::FATAL,
           "Cannot allocate memory for driver");

    // Solve
    driver->execute();

  } // Catch and handle Quina::Exceptions
    catch (Exception& qe) {
      error = qe.handleException(driver);
    }
    // Catch std::exceptions and transform them into Quinoa::Exceptions without
    // file:line:func information
    catch (std::exception& se) {
      Exception qe(ExceptType::RUNTIME, se.what());
      error = qe.handleException(driver);
    }
    // Catch uncaught exceptions and still do cleanup
    catch (...) {
      Exception qe(ExceptType::UNCAUGHT, "Non-standard exception");
      error = qe.handleException(driver);
    }

  // Finalize and deallocate driver and memory manager
  if (driver) { delete driver; driver = nullptr; }
  if (memory) { delete memory; memory = nullptr; }

  // Return error code
  return static_cast<int>(error);
}
