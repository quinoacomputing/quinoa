//******************************************************************************
/*!
  \file      src/Mesh/MeshException.C
  \author    J. Bakosi
  \date      Fri Apr 26 12:57:11 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshException class definition
  \details   MeshException class definition
*/
//******************************************************************************

#include <MeshException.h>

using namespace Quinoa;

ErrCode
MeshException::handleException(Driver* const driver)
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
