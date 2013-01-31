//******************************************************************************
/*!
  \file      src/Parser/ParserException.C
  \author    J. Bakosi
  \date      Wed 30 Jan 2013 07:02:52 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ParserException
  \details   ParserException
*/
//******************************************************************************

#include <ParserException.h>

using namespace Quinoa;

ErrCode
ParserException::handleException(Driver* driver)
//******************************************************************************
//  Handle ParserException
//! \author J. Bakosi
//******************************************************************************
{
  // Start error message
  m_message = ParserMsg[static_cast<int>(m_except)];
  if (m_throwerMsg.size()) m_message += m_throwerMsg;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
