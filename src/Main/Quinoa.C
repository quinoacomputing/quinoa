//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Thu Aug 29 15:21:58 2013
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

using namespace quinoa;

namespace quinoa {

static void echoName()
//******************************************************************************
//  Echo Name
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << "=========================================\n";
  std::cout << "Quinoa: Lagrangian particle hydrodynamics\n";
  std::cout << "=========================================" << std::endl;
}

static void echoBuildInfo()
//******************************************************************************
//  Echo build environment
//! \details Echo information read from [build]/Base/QuinoaConfig.h filled by
//!          CMake based on src/MainQuinoaConfig.h.in
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << "\nBuild environment:"
               "\n------------------\n";
  std::cout << " * Executable                  : " << QUINOA_EXECUTABLE << "\n";
  std::cout << " * Version                     : " << QUINOA_VERSION << "\n";
  std::cout << " * Release                     : " << QUINOA_RELEASE << "\n";
  std::cout << " * Revision                    : " << QUINOA_GIT_COMMIT << "\n";
  std::cout << " * CMake build type            : " << QUINOA_BUILD_TYPE << "\n";
  std::cout << " * MPI C++ compiler            : " << QUINOA_MPI_COMPILER<<"\n";
  std::cout << " * MPI underlying C++ compiler : " << QUINOA_COMPILER << "\n";
  std::cout << " * Build date                  : " << QUINOA_BUILD_DATE << "\n";
#ifdef NDEBUG
  std::cout << " * Built without asserts" << "\n";
#else  // NDEBUG
  std::cout << " * Built with asserts" << "\n";
#endif // NDEBUG
  std::cout << std::endl;
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

    // Echo program name
    echoName();
    // Echo build environment
    echoBuildInfo();

    // Query, setup, and echo parallel enviroment
    Paradigm paradigm;
    paradigm.echo();

    // Initialize memory manager
    memory = new (std::nothrow) Memory(&paradigm);
    ErrChk(memory != nullptr, ExceptType::FATAL,
           "No memory for a memory manager?");

    // Allocate and initialize driver
    driver = new (std::nothrow) QuinoaDriver(argc, argv, memory, &paradigm);
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
