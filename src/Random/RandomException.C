//******************************************************************************
/*!
  \file      src/Random/RandomException.C
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 04:18:40 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RandomException class definition
  \details   RandomException class definition
*/
//******************************************************************************

#include <RandomException.h>

using namespace Quinoa;

ErrCode
RandomException::handleException(Driver* driver)
//******************************************************************************
//  Handle RandomException
//! \author J. Bakosi
//******************************************************************************
{
  // Contribute to error message
  message = RndMsg[static_cast<int>(m_except)] + message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
