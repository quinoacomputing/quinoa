//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Thu 27 Jun 2013 09:20:45 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>

#include <QuinoaConfig.h>
#include <Paradigm.h>
#include <Memory.h>
#include <Driver.h>

using namespace std;
using namespace Quinoa;

namespace Quinoa {

static void echoName()
//******************************************************************************
//  Echo Name
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "=========================================\n";
  cout << "Quinoa: Lagrangian particle hydrodynamics\n";
  cout << "=========================================" << endl;
}

static void echoBuildInfo()
//******************************************************************************
//  Echo build environment
//! \details Echo information read from [build]/Base/QuinoaConfig.h filled by
//!          CMake based on src/MainQuinoaConfig.h.in
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "\nBuild environment:"
          "\n------------------\n";
  cout << " * Executable                  : " << QUINOA_EXECUTABLE << "\n";
  cout << " * Version                     : " << QUINOA_VERSION << "\n";
  cout << " * Release                     : " << QUINOA_RELEASE << "\n";
  cout << " * Revision                    : " << QUINOA_GIT_COMMIT << "\n";
  cout << " * CMake build type            : " << QUINOA_BUILD_TYPE << "\n";
  cout << " * MPI C++ compiler            : " << QUINOA_MPI_COMPILER << "\n";
  cout << " * MPI underlying C++ compiler : " << QUINOA_COMPILER << "\n";
  cout << " * Build date                  : " << QUINOA_BUILD_DATE << "\n";
#ifdef NDEBUG
  cout << " * Built without asserts" << "\n";
#else  // NDEBUG
  cout << " * Built with asserts" << "\n";
#endif // NDEBUG
  cout << endl;
}

} // namespace Quinoa

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa main
//! \details Quinoa main -- This is where everything starts...
//! \author  J. Bakosi
//******************************************************************************
{
  Memory* memory = nullptr;
  Driver* driver = nullptr;

  ErrCode error = ErrCode::HAPPY;
  try {

    // Echo program name
    echoName();
    // Echo build environment
    echoBuildInfo();

    // Query, setup, and echo parallel enviroment
    Paradigm paradigm;
    paradigm.echo();

    // Initialize memory manager
    memory = new (nothrow) Memory(&paradigm);
    ErrChk(memory != nullptr, ExceptType::FATAL,
           "No memory for a memory manager?");

    // Allocate and initialize driver
    driver = new (nothrow) Driver(argc, argv, memory, &paradigm);
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
    catch (exception& se) {
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
