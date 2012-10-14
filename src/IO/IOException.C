//******************************************************************************
/*!
  \file      src/IO/IOException.C
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 05:26:44 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     IOException class definition
  \details   IOException class definition
*/
//******************************************************************************

#include <IOException.h>

using namespace Quinoa;

ErrCode
IOException::handleException(Driver* driver)
//******************************************************************************
//  Handle IOException
//! \author J. Bakosi
//******************************************************************************
{
  // Start error message
  message = IOMsg[static_cast<Int>(m_except)];
  if (m_filename.size()) message += m_filename;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
