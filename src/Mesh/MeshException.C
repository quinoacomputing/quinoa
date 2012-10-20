//******************************************************************************
/*!
  \file      src/Mesh/MeshException.C
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 04:18:26 PM MDT
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
  message = MeshMsg[static_cast<int>(m_except)];
  if (m_throwerMsg.size()) message += m_throwerMsg;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
