//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 05:41:39 AM MST
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
      cerr << "WARNING: " << m_message << endl;
      return NONFATAL;
    case ERROR:
      cerr << "ERROR: " << m_message << endl;
      return NONFATAL;
    case UNCAUGHT:  // Warn and fall through FATAL
      cerr << "UNCAUGHT EXCEPTION: " << m_message << endl;
    case FATAL:     // Attempt cleanup and exit
      cerr << "FATAL ERROR: " << m_message << endl
           << "Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;
    default:
      return NONFATAL;
  }
}
