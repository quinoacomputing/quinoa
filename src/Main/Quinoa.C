//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Mon Oct  7 10:58:37 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <Base.h>
#include <Handler.h>
#include <Paradigm.h>
#include <QuinoaDriver.h>
#include <QuinoaPrint.h>

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

    // Create the essentials
    InputDeck control;                       //!< Control
    QuinoaPrint print(control);              //!< Pretty printer
    Paradigm paradigm(print);                //!< Parallel compute environment
    Timer timer;                             //!< Timer

    // Bundle up essentials
    Base base(print, paradigm, control, timer);

    // Echo program name
    echoHeader(print, "Quinoa: Lagrangian particle hydrodynamics");

    // Echo environment
    print.part("Environment");
    echoBuildEnv(print, QUINOA_EXECUTABLE);  //!< Build environment
    paradigm.echo();                         //!< Parallel compute env
    echoRunEnv(print, argc, argv);           //!< Runtime environment

    // Create driver
    QuinoaDriver driver(argc, argv, base);

    // Execute
    driver.execute();

  } catch (...) {
      error = processException();
    }

  // Return error code
  return static_cast<int>(error);
}
