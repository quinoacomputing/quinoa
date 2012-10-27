//******************************************************************************
/*!
  \file      src/Statistics/StatException.C
  \author    J. Bakosi
  \date      Sat 27 Oct 2012 11:38:38 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistics exception
  \details   Statistics exception
*/
//******************************************************************************

#include <StatException.h>

using namespace Quinoa;

ErrCode
StatException::handleException(Driver* driver)
//******************************************************************************
//  Handle StatException
//! \author J. Bakosi
//******************************************************************************
{
  // Start error message
  message = StatMsg[static_cast<int>(m_except)];

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
