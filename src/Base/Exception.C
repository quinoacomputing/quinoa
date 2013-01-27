//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Sun 27 Jan 2013 11:37:50 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class definition
  \details   Exception base class definition
*/
//******************************************************************************

#include <iostream>
#include <sstream>

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
  // Add "file:line:func" to message
  stringstream s;
  if (m_except != UNCAUGHT) {
    s << "<< " << m_message << " >> Exception in " << m_file << ":" << m_line
      << ": " << m_func << "\n";
  }
  m_message = s.str();

  switch (m_except) {

    case WARNING:
      cout << "WARNING: " << m_message;
      return NONFATAL;

    case ERROR:
      cout << "ERROR: " << m_message;
      return NONFATAL;

    case UNCAUGHT:  // Warn and return FATAL
      cout << "UNKNOWN EXCEPTION\n"
           << "Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;

    case FATAL:     // Attempt cleanup and exit
      cout << "FATAL ERROR: " << m_message
           << "Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;

    default:
      return NONFATAL;
  }
}
