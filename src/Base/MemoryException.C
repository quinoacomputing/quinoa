//******************************************************************************
/*!
  \file      src/Base/MemoryException.C
  \author    J. Bakosi
  \date      Fri Sep 14 17:32:49 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MemoryException class definition
  \details   MemoryException class definition
*/
//******************************************************************************

#include <iostream>

#include <MemoryException.h>

using namespace Quinoa;

ErrCode
MemoryException::handleException(Driver* driver)
//******************************************************************************
//  Handle MemoryException
//! \author J. Bakosi
//******************************************************************************
{
  // Output message
  cerr << "Memory exception: " << MemMsg[static_cast<Int>(m_except)] << endl;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
