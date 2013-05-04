//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Sat 04 May 2013 06:56:21 AM MDT
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
                     const std::string& message,
                     int number) noexcept
//******************************************************************************
//  Constructor: generate error message
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
try :
  m_except(move(except)),
  m_message(move(message))
{

  stringstream s;
  s << m_message;
  if (number) s << number;
  m_message = s.str();

} // Catch std::exception
  catch (exception& se) {
    // Emit warning and continue
    cout << "RUNTIME ERROR in Exception constructor: " << se.what() << endl
         << "Continuing anyway..." << endl;
  }
  // Catch uncaught exceptions
  catch (...) {
    // Emit warning and continue
    cout << "UNKNOWN EXCEPTION in Exception constructor" << endl
         << "Continuing anyway..." << endl;
  }

void
Exception::augment(const std::string& message) noexcept
//******************************************************************************
//  Augment message after its construction
//! \details ICC: in C++11 after initializer lists are properly supported,
//! VSLException's constructor can be made empty, invalidating this function.
//! No-throw guarantee: this member function never throws exceptions.
//! \param[in]  message  Message to add to end of message
//! \author J. Bakosi
//******************************************************************************
{
  try {

    m_message += message;

    } // Catch std::exception
      catch (exception& se) {
        // Emit warning and continue
        cout << "RUNTIME ERROR in Exception::augment(): " << se.what() << endl
             << "Continuing anyway..." << endl;
      }
      // Catch uncaught exceptions
      catch (...) {
        // Emit warning and continue
        cout << "UNKNOWN EXCEPTION in Exception::augment()" << endl
             << "Continuing anyway..." << endl;
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
