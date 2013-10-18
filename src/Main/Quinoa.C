//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Fri Oct 18 11:36:56 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <Handler.h>
#include <Print.h>
#include <QuinoaDriver.h>

using namespace quinoa;
using namespace tk;

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Lagrangian particle hydrodynamics
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
    echoHeader(print, "Quinoa: Lagrangian particle hydrodynamics");

    // Echo environment
    print.part("Environment");
    echoBuildEnv(print, QUINOA_EXECUTABLE);  //!< Build environment
    echoRunEnv(print, argc, argv);           //!< Runtime environment

    // Create driver
    QuinoaDriver driver(argc, argv, print);

    // Execute
    driver.execute();

  } catch (...) {
      error = processException();
    }

  // Return error code
  return static_cast<int>(error);
}
