//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Fri Oct 18 11:38:13 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa random number generator test suite
  \details   Quinoa random number generator test suite
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <Handler.h>
#include <Print.h>
#include <RNGTestDriver.h>

using namespace rngtest;
using namespace tk;

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Random number generator test suite
//! \author  J. Bakosi
//******************************************************************************
{
  ErrCode error = ErrCode::SUCCESS;

  try {

    // Install our own new-handler
    std::set_new_handler(newHandler);
    // Install our own terminate-handler
    std::set_terminate(terminateHandler);
    // Install our own unexpected-handler
    std::set_unexpected(unexpectedHandler);

    // Create pretty printer
    Print print;

    // Echo program name
    echoHeader(print, "Quinoa: Random number generator test suite");

    // Echo environment
    print.part("Environment");
    echoBuildEnv(print, RNGTEST_EXECUTABLE);  //!< Build environment
    echoRunEnv(print, argc, argv);            //!< Runtime environment

    // Create driver
    RNGTestDriver driver(argc, argv, print);

    // Execute
    driver.execute();

  } catch (...) {
      error = processException();
    }

  // Return error code
  return static_cast<int>(error);
}
