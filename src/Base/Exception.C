//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 07:25:32 PM MDT
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
    case WARNING:
      cerr << "WARNING: " << message << endl;
      return NONFATAL;
    case ERROR:
      cerr << "ERROR: " << message << endl;
      return NONFATAL;
    case UNCAUGHT:  // Warn and fall through FATAL
      cerr << "UNCAUGHT EXCEPTION: " << message << endl;
    case FATAL:     // Attempt cleanup and exit
      cerr << "FATAL ERROR: " << message << endl
           << "Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;
    default:
      return NONFATAL;
  }
}
