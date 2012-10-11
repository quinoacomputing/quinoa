//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Wed 10 Oct 2012 01:42:29 PM EDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class definition
  \details   Exception base class definition
*/
//******************************************************************************

#include <iostream>

#include <Driver.h>
#include <Exception.h>

using namespace Quinoa;

ErrCode
Exception::handleException(Driver* driver)
//******************************************************************************
//  Handle Exception (criticality)
//! \author J. Bakosi
//******************************************************************************
{
  switch (m_except) {
    // If exception is uncaught, warn and fall through FATAL
    case UNCAUGHT:
      cerr << "WARNING: Uncaught exception" << endl;
    // If exception is fatal, do cleanup and exit
    case FATAL:
      cerr << "FATAL ERROR: Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;
    default:
      return NONFATAL;
  }
}
