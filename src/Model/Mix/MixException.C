//******************************************************************************
/*!
  \file      src/Model/Mix/MixException.C
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 11:26:25 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model exception handler
  \details   Mix model exception handler
*/
//******************************************************************************

#include <MixException.h>

using namespace Quinoa;

ErrCode
MixException::handleException(Driver* driver)
//******************************************************************************
//  Handle MixException
//! \author J. Bakosi
//******************************************************************************
{
  // Contribute to error message
  m_message = MixMsg[static_cast<int>(m_except)] + m_message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
