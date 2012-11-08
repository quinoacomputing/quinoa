//******************************************************************************
/*!
  \file      src/Base/MemoryException.C
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 05:20:40 AM MST
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
  m_message = MemMsg[static_cast<int>(m_except)];

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
