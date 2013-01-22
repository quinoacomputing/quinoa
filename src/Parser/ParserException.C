//******************************************************************************
/*!
  \file      src/Parser/ParserException.C
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 08:37:49 PM MST
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
  // Contribute to error message
  m_message = ParserMsg[static_cast<int>(m_except)] + m_message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
