//******************************************************************************
/*!
  \file      src/Random/MKLException.C
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 05:28:18 PM MDT
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
  // Start error message
  message = MKLMsg[static_cast<Int>(m_except)];

  // Handle Exception (criticality)
  return RandomException::handleException(driver);
}
