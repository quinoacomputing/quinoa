//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Mon Apr 29 16:00:20 2013
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

Exception::Exception(const ExceptType except,
                     const string& message,
                     const string& file,
                     const string& func,
                     unsigned int line) noexcept :
  m_file(file),
  m_func(func),
  m_line(line),
  m_message(message),
  m_except(except)
//******************************************************************************
//  Constructor: generate error message
//! \author J. Bakosi
//******************************************************************************
{
  if (m_line) {         // append information in file:line:func if available
    stringstream s;
    s << m_message << endl
      << ">>> Exception in " << m_file << ":" << m_line << ": " << m_func;
    m_message = s.str();
  }
}

ErrCode
Exception::handleException(Driver* const driver)
//******************************************************************************
//  Handle Exception: Print cumulative message and handle criticality
//! \author J. Bakosi
//******************************************************************************
{
  switch (m_except) {

    case WARNING:
      cout << ">>> WARNING: " << what() << endl;
      return NONFATAL;

    case CUMULATIVE:
      cout << ">>> CUMULATIVE ERROR: " << what() << endl;
      return NONFATAL;

    case ERROR:
      cout << ">>> ERROR: " << what() << endl;
      return NONFATAL;

    case FATAL:
      cout << ">>> FATAL ERROR: " << what() << endl
           << ">>> Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;

    case RUNTIME:
      cout << ">>> RUNTIME ERROR: " << what() << endl
           << ">>> Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;

    case UNCAUGHT:
      cout << ">>> UNKNOWN ERROR\n"
           << ">>> Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return FATAL_ERROR;

    default:
      return NONFATAL;
  }
}
