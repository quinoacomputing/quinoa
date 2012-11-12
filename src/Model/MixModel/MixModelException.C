//******************************************************************************
/*!
  \file      src/Model/MixModel/MixModelException.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 09:34:14 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model exception handler
  \details   Mix model exception handler
*/
//******************************************************************************

#include <MixModelException.h>

using namespace Quinoa;

ErrCode
MixModelException::handleException(Driver* driver)
//******************************************************************************
//  Handle MixModelException
//! \author J. Bakosi
//******************************************************************************
{
  // Contribute to error message
  m_message = MixModelMsg[static_cast<int>(m_except)] + m_message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
