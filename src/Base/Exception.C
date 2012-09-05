//******************************************************************************
/*!
  \file      src/Base/Exception.C
  \author    J. Bakosi
  \date      Tue 04 Sep 2012 10:45:42 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Exception base class definition
  \details   Exception base class definition
*/
//******************************************************************************

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
  // If exception is fatal, do cleanup and exit
  if (m_exception == FATAL) {
    driver->finalize();
    return FATAL_ERROR;
  } else {
    return NONFATAL_ERROR;
  }
}
