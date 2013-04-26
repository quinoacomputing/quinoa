//******************************************************************************
/*!
  \file      src/Physics/PhysicsException.C
  \author    J. Bakosi
  \date      Fri Apr 26 12:57:51 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     PhysicsException
  \details   PhysicsException
*/
//******************************************************************************

#include <PhysicsException.h>

using namespace Quinoa;

ErrCode
PhysicsException::handleException(Driver* const driver)
//******************************************************************************
//  Handle PhysicsException
//! \author J. Bakosi
//******************************************************************************
{
  // Contribute to error message
  m_message = PhysicsMsg[static_cast<int>(m_except)] + m_message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
