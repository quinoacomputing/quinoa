//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 12:20:27 PM MST
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

static void echoBuildInfo()
//******************************************************************************
//  Echo build environment
//! \details Echo information read from \<build\>/Base/QuinoaConfig.h filled by
//!          CMake based on \<source\>/MainQuinoaConfig.h.in
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "\nBuild environment:\n";
  cout << " * Executable         : " << QUINOA_EXECUTABLE << "\n";
  cout << " * Version            : " << QUINOA_VERSION << "\n";
  cout << " * Release            : " << QUINOA_RELEASE << "\n";
  cout << " * Git commit         : " << QUINOA_GIT_COMMIT << "\n";
  cout << " * Configuration      : " << QUINOA_CONFIGURATION << "\n";
  cout << " * Third-party prefix : " << QUINOA_THIRD_PARTY_PREFIX << "\n";
  cout << " * Compiler           : " << QUINOA_COMPILER << "\n";
  cout << " * Build type         : " << QUINOA_BUILD_TYPE;
#ifdef NDEBUG
  cout << " (no error checking, no exception handling)" << "\n";
#else  // NDEBUG
  cout << " (with error checking and exception handling)" << "\n";
#endif // NDEBUG
  cout << " * Build date         : " << QUINOA_BUILD_DATE << "\n";
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
  // Echo build environment
  echoBuildInfo();

  // Query and echo parallel enviroment
  Paradigm paradigm;
  paradigm.echo();

  // Get number of OpenMP threads
  const OpenMP* omp = paradigm.getOpenMP();
  const int ntomp = omp->nthread();

  // Initialize memory manager and driver
  Memory memStore(ntomp);
  Driver driver(&memStore);

  ErrCode error = NO_ERROR;
  try {

    driver.setup();
    driver.solve();

  } catch (Exception& e) { error = e.handleException(&driver); }
    catch (...) { // catch uncaught exceptions
      Exception e(UNCAUGHT);
      error = e.handleException(&driver);
    }

  //!< Finalize
  driver.finalize();

  if (error != FATAL_ERROR) cout << "Normal finish." << endl;
}
