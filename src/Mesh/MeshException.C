//******************************************************************************
/*!
  \file      src/Mesh/MeshException.C
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 07:20:41 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshException class definition
  \details   MeshException class definition
*/
//******************************************************************************

#include <MeshException.h>

using namespace Quinoa;

ErrCode
MeshException::handleException(Driver* driver)
//******************************************************************************
//  Handle MeshException
//! \author J. Bakosi
//******************************************************************************
{
  // Start error message
  message = MeshMsg[static_cast<Int>(m_except)];
  if (m_throwerMsg.size()) message += m_throwerMsg;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
