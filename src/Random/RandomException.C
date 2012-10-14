//******************************************************************************
/*!
  \file      src/Random/RandomException.C
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 05:28:33 PM MDT
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
  message = RndMsg[static_cast<Int>(m_except)] + message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
