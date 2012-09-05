//******************************************************************************
/*!
  \file      src/Base/MemoryException.C
  \author    J. Bakosi
  \date      Tue 04 Sep 2012 10:48:48 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MemoryException class definition
  \details   MemoryException class definition
*/
//******************************************************************************

#include <iostream>

#include <MemoryException.h>

using namespace Quinoa;

ErrorCode
MemoryException::handleException(Driver* driver)
//******************************************************************************
//  Handle MemoryException
//! \author J. Bakosi
//******************************************************************************
{
  // Output message
  cerr << "Memory exception: " << MemoryMessage[m_exception] << endl;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
