//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Sun 27 Jan 2013 12:57:10 PM MST
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

string
Exception::genericWhat()
//******************************************************************************
//  Generate generic exception message
//! \author J. Bakosi
//******************************************************************************
{
  // Add "file:line:func" to cumulative message
  stringstream s;
  s << m_message << ", Exception in "
    << m_file << ":" << m_line << ": " << m_func << "\n";
  return s.str();
}

string
Exception::runtimeWhat()
//******************************************************************************
//  Generate std::runtime_error exception message
//! \author J. Bakosi
//******************************************************************************
{
  // "file:line:func" uninitialized, but have a message
  stringstream s;
  s << m_message << "\n";
  return s.str();
}

ErrCode
Exception::handleException(Driver* driver)
//******************************************************************************
//  Handle Exception: Print cumulative message and handle criticality
//! \author J. Bakosi
//******************************************************************************
{
  switch (m_except) {

    case WARNING:
      cout << "WARNING: " << genericWhat();
      return NONFATAL;

    case CUMULATIVE:
      cout << "CUMULATIVE ERROR: " << genericWhat();
      return NONFATAL;

    case ERROR:
      cout << "ERROR: " << genericWhat();
      return NONFATAL;

    case FATAL:
      cout << "FATAL ERROR: " << genericWhat()
           << "Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;

    case RUNTIME:
      cout << "RUNTIME ERROR: " << runtimeWhat()
           << "Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;

    case UNCAUGHT:
      cout << "UNKNOWN ERROR\n"
           << "Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;

    default:
      return NONFATAL;
  }
}
