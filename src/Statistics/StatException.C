//******************************************************************************
/*!
  \file      src/Statistics/StatException.C
  \author    J. Bakosi
  \date      Sat 16 Mar 2013 09:44:17 AM MDT
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
  if (m_throwerMsg.size()) m_message += m_throwerMsg;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
