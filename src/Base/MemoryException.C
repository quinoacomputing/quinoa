//******************************************************************************
/*!
  \file      src/Base/MemoryException.C
  \author    J. Bakosi
  \date      Fri Apr 26 12:56:28 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MemoryException class definition
  \details   MemoryException class definition
*/
//******************************************************************************

#include <MemoryException.h>

using namespace Quinoa;

ErrCode
MemoryException::handleException(Driver* const driver)
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
