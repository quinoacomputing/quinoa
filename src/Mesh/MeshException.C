//******************************************************************************
/*!
  \file      src/Mesh/MeshException.C
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 05:21:01 AM MST
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
  m_message = MeshMsg[static_cast<int>(m_except)];
  if (m_throwerMsg.size()) m_message += m_throwerMsg;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
