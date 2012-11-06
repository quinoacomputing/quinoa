//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Mon 05 Nov 2012 06:42:12 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <iostream>

#include <QuinoaConfig.h>
#include <Paradigm.h>
#include <Setup.h>

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

}
