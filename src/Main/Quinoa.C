//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Mon Jan 28 08:03:33 2013
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
  cout << "===============================================\n";
  cout << "Quinoa: Lagrangian particle hydrodynamics\n";
  cout << "===============================================" << endl;
}

static void echoBuildInfo()
//******************************************************************************
//  Echo build environment
//! \details Echo information read from \<build\>/Base/QuinoaConfig.h filled by
//!          CMake based on \<source\>/MainQuinoaConfig.h.in
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "\nBuild environment:\n";
  cout << " * Executable        : " << QUINOA_EXECUTABLE << "\n";
  cout << " * Version           : " << QUINOA_VERSION << "\n";
  cout << " * Release           : " << QUINOA_RELEASE << "\n";
  cout << " * Revision          : " << QUINOA_GIT_COMMIT << "\n";
  cout << " * Configuration     : " << QUINOA_CONFIGURATION << "\n";
  cout << " * Third-party prefix: " << QUINOA_THIRD_PARTY_PREFIX << "\n";
  cout << " * Compiler          : " << QUINOA_COMPILER << "\n";
  cout << " * Build type        : " << QUINOA_BUILD_TYPE;
#ifdef NDEBUG
  cout << " (no error checking, no exception handling)" << "\n";
#else  // NDEBUG
  cout << " (with error checking and exception handling)" << "\n";
#endif // NDEBUG
  cout << " * Build date        : " << QUINOA_BUILD_DATE << "\n";
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
  // Echo name
  echoName();

  // Echo build environment
  echoBuildInfo();

  // Query and echo parallel enviroment
  Paradigm paradigm;
  paradigm.echo();

  // Initialize memory manager and driver
  Memory memory(&paradigm);
  Driver driver(argc, argv, &memory, &paradigm);

  ErrCode error = HAPPY;
  // This main try-catch block is executed even in NDEBUG mode, otherwise
  // memory may leak
  try {

    driver.setup();
    driver.solve();

  } // Catch and handle Quina::Exceptions
    catch (Exception& qe) {
      error = qe.handleException(&driver);
    }
    // Catch std::exceptions and transform them into Quinoa::Exceptions
    catch (exception& se) {
      Exception qe(RUNTIME, se.what());
      error = qe.handleException(&driver);
    }
    // Catch uncaught exceptions and still do cleanup
    catch (...) {
      Exception qe(UNCAUGHT);
      error = qe.handleException(&driver);
    }

  // Finalize
  driver.finalize();

  // Return the error code
  return error;
}
