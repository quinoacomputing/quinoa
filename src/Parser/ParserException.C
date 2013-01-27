//******************************************************************************
/*!
  \file      src/Parser/ParserException.C
  \author    J. Bakosi
  \date      Sun 27 Jan 2013 11:10:06 AM MST
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

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
