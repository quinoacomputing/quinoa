//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Fri 09 Nov 2012 06:56:37 PM MST
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

int main(int argc, char* argv[]) {

  // Echo build environment
  cout << endl << "Build environment:" << endl;
  cout << " * Executable         : " << QUINOA_EXECUTABLE << endl;
  cout << " * Version            : " << QUINOA_VERSION << endl;
  cout << " * Release            : " << QUINOA_RELEASE << endl;
  cout << " * Git commit         : " << QUINOA_GIT_COMMIT << endl;
  cout << " * Configuration      : " << QUINOA_CONFIGURATION << endl;
  cout << " * Third-party prefix : " << QUINOA_THIRD_PARTY_PREFIX << endl;
  cout << " * Compiler           : " << QUINOA_COMPILER << endl;
  cout << " * Build type         : " << QUINOA_BUILD_TYPE << endl;
  cout << " * Build date         : " << QUINOA_BUILD_DATE << endl;
  cout << endl;

  // Query and echo parallel programming enviroment
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
