//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 07:48:38 PM MST
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
  s << "<< " << m_message << " >> Exception in " << m_file << ":" << m_line
    << ": " << m_func << "\n";
  m_message = s.str();

  switch (m_except) {
    case WARNING:
      cout << "WARNING: " << m_message;
      return NONFATAL;
    case ERROR:
      cout << "ERROR: " << m_message;
      return NONFATAL;
    case UNCAUGHT:  // Warn and fall through FATAL
      cout << "UNCAUGHT EXCEPTION" << m_message;
    case FATAL:     // Attempt cleanup and exit
      cout << "FATAL ERROR: " << m_message
           << "Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;
    default:
      return NONFATAL;
  }
}
