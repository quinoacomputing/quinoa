//******************************************************************************
/*!
  \file      src/IO/IOException.C
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 04:12:34 PM MDT
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
  message = IOMsg[static_cast<int>(m_except)];
  if (m_filename.size()) message += m_filename;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
