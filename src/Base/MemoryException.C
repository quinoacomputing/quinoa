//******************************************************************************
/*!
  \file      src/Base/MemoryException.C
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 04:11:26 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MemoryException class definition
  \details   MemoryException class definition
*/
//******************************************************************************

#include <MemoryException.h>

using namespace Quinoa;

ErrCode
MemoryException::handleException(Driver* driver)
//******************************************************************************
//  Handle MemoryException
//! \author J. Bakosi
//******************************************************************************
{
  // Start error message
  message = MemMsg[static_cast<int>(m_except)];

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
