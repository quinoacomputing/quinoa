//******************************************************************************
/*!
  \file      src/IO/IOException.C
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 05:20:50 AM MST
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
  m_message = IOMsg[static_cast<int>(m_except)];
  if (m_filename.size()) m_message += m_filename;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
