//******************************************************************************
/*!
  \file      src/Random/RandomException.C
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 05:21:11 AM MST
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
  m_message = RndMsg[static_cast<int>(m_except)] + m_message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
