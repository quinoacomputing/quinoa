//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Mon May  6 15:20:44 2013
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
//! \details Echo information read from \<build\>/Base/QuinoaConfig.h filled by
//!          CMake based on \<source\>/MainQuinoaConfig.h.in
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "\nBuild environment:"
          "\n------------------\n";
  cout << " * Executable        : " << QUINOA_EXECUTABLE << "\n";
  cout << " * Version           : " << QUINOA_VERSION << "\n";
  cout << " * Release           : " << QUINOA_RELEASE << "\n";
  cout << " * Revision          : " << QUINOA_GIT_COMMIT << "\n";
  cout << " * Configuration     : " << QUINOA_CONFIGURATION << "\n";
  cout << " * Third-party prefix: " << QUINOA_THIRD_PARTY_PREFIX << "\n";
  cout << " * Compiler          : " << QUINOA_COMPILER << "\n";
  cout << " * Build type        : " << QUINOA_BUILD_TYPE;
#ifdef NDEBUG
  cout << " (no asserts)" << "\n";
#else  // NDEBUG
  cout << " (with asserts)" << "\n";
#endif // NDEBUG
  cout << " * Build date        : " << QUINOA_BUILD_DATE << "\n";
  cout << endl;
}

} // namespace Quinoa

static void finalize(const Driver* driver, const Memory* memory) noexcept
//******************************************************************************
//  Finalize for main()
//! \details This is to avoid code-duplication in main() due to exception
//!          handling.
//!          Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
 if (driver) { delete driver; driver = nullptr; }
 if (memory) { delete memory; memory = nullptr; }
}

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa main
//! \details Quinoa main -- This is where everything starts...
//! \author  J. Bakosi
//******************************************************************************
{
  Driver* driver = nullptr;
  Memory* memory = nullptr;

  ErrCode error = HAPPY;
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
    Errchk(memory != nullptr, FATAL, "No memory for a memory manager?");

    // Allocate and initialize driver
    driver = new (nothrow) Driver(argc, argv, memory, &paradigm);
    Errchk(driver != nullptr, FATAL, "Cannot allocate memory for driver");

    // Setup and solve
    driver->setup();
    driver->solve();

  } // Catch and handle Quina::Exceptions
    catch (Exception& qe) {
      finalize(driver, memory);
      error = qe.handleException(driver);
    }
    // Catch std::exceptions and transform them into Quinoa::Exceptions without
    // file:line:func information
    catch (exception& se) {
      finalize(driver, memory);
      Exception qe(RUNTIME, se.what());
      error = qe.handleException(driver);
    }
    // Catch uncaught exceptions and still do cleanup
    catch (...) {
      finalize(driver, memory);
      Exception qe(UNCAUGHT);
      error = qe.handleException(driver);
    }

  // Finalize and deallocate driver and memory manager
  finalize(driver, memory);

  // Return error code
  return error;
}
