//******************************************************************************
/*!
  \file      src/Random/MKLException.C
  \author    J. Bakosi
  \date      Wed Nov  7 17:42:14 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKLException class definition
  \details   MKLException class definition
*/
//******************************************************************************

#include <MKLException.h>

using namespace Quinoa;

ErrCode
MKLException::handleException(Driver* driver)
//******************************************************************************
//  Handle MKLException
//! \author J. Bakosi
//******************************************************************************
{
  // Contribute to error message
  m_message = MKLMsg[static_cast<int>(m_except)] + m_message;

  // Handle Exception (criticality)
  return RandomException::handleException(driver);
}
