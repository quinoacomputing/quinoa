//******************************************************************************
/*!
  \file      src/Statistics/StatException.C
  \author    J. Bakosi
  \date      Fri Apr 26 12:58:25 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistics exception
  \details   Statistics exception
*/
//******************************************************************************

#include <StatException.h>

using namespace Quinoa;

ErrCode
StatException::handleException(Driver* const driver)
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
