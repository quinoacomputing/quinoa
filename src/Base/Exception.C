//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Fri Sep 14 17:27:05 2012
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
    case ExceptType::UNCAUGHT:
      cerr << "WARNING: Uncaught exception" << endl;
    // If exception is fatal, do cleanup and exit
    case ExceptType::FATAL:
      cerr << "FATAL ERROR: Attempting cleanup & graceful exit..." << endl;
      driver->finalize();
      return ErrCode::FATAL;
    default:
      return ErrCode::NONFATAL;
  }
}
