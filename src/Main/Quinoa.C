//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Sun 04 Nov 2012 10:05:54 PM MST
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
  cout << "Build environment:" << endl;
  cout << " * Executable         : " << QUINOA_EXECUTABLE << endl;
  cout << " * Version            : " << QUINOA_VERSION << endl;
  cout << " * Release            : " << QUINOA_RELEASE << endl;
  cout << " * Revision           : " << QUINOA_REVISION << endl;
  cout << " * Revision date      : " << QUINOA_REVISION_DATE << endl;
  cout << " * Configuration      : " << QUINOA_CONFIGURATION << endl;
  cout << " * Third-party prefix : " << QUINOA_THIRD_PARTY_PREFIX << endl;
  cout << " * Compiler           : " << QUINOA_COMPILER << endl;
  cout << " * Build type         : " << QUINOA_BUILD_TYPE << endl;
  cout << " * Build date         : " << QUINOA_BUILD_DATE << endl;
  cout << endl;

  // Query parallel programming enviroment
  Paradigm paradigm;
  paradigm.echo();

}
