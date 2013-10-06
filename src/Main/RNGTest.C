//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 02:36:40 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa random number generator test suite
  \details   Quinoa random number generator test suite
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <Base.h>
#include <Handler.h>
#include <Paradigm.h>
#include <RNGTestDriver.h>
#include <RNGTestPrint.h>

using namespace quinoa;
using namespace rngtest;

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

    // Create the essentials
    rngtest::InputDeck control;             //!< Control
    RNGTestPrint print(control);            //!< Pretty printer
    Paradigm paradigm(print);               //!< Parallel compute environment
    Timer timer;                            //!< Timer

    // Bundle up essentials
    rngtest::Base base(print, paradigm, control, timer);

    // Echo program name
    init::echoHeader(print, "Quinoa: Random number generator test suite");

    // Echo environment
    print.part("Environment");
    init::echoBuildEnv(print, RNGTEST_EXECUTABLE);  //!< Build environment
    paradigm.echo();                                //!< Parallel compute env
    init::echoRunEnv(print, argc, argv);            //!< Runtime environment

    // Create driver
    RNGTestDriver driver(argc, argv, base);

    // Execute
    driver.execute();

  } catch (...) {
      error = processException();
    }

  // Return error code
  return static_cast<int>(error);
}
