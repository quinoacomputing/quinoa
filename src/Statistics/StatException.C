//******************************************************************************
/*!
  \file      src/Statistics/StatException.C
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 05:21:49 AM MST
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
  m_message = StatMsg[static_cast<int>(m_except)];

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
