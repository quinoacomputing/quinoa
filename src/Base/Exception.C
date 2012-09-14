//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Fri Sep 14 15:56:38 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class definition
  \details   Exception base class definition
*/
//******************************************************************************

#include <iostream>

#include <Driver.h>
#include <Exception.h>

using namespace Quinoa;

ErrorCode
Exception::handleException(Driver* driver)
//******************************************************************************
//  Handle Exception (criticality)
//! \author J. Bakosi
//******************************************************************************
{
  switch (m_exception) {
    case UNCAUGHT: // If exception is uncaught, warn and fall through FATAL
      cerr << "WARNING: Uncaught exception" << endl;
    case FATAL:    // If exception is fatal, do cleanup and exit
      cerr << "FATAL ERROR: Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;
    default:
      return NONFATAL_ERROR;
  }
}
