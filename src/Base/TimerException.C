//******************************************************************************
/*!
  \file      src/Base/TimerException.C
  \author    J. Bakosi
  \date      Fri Apr 26 12:58:51 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TimerException class definition
  \details   TimerException class definition
*/
//******************************************************************************

#include <TimerException.h>

using namespace Quinoa;

ErrCode
TimerException::handleException(Driver* const driver)
//******************************************************************************
//  Handle TimerException
//! \author J. Bakosi
//******************************************************************************
{
  // Start error message
  m_message = TimerMsg[static_cast<int>(m_except)];

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
