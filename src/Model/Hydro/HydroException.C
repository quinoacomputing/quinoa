//******************************************************************************
/*!
  \file      src/Model/Hydro/HydroException.C
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 11:31:54 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model exception handler
  \details   Hydro model exception handler
*/
//******************************************************************************

#include <HydroException.h>

using namespace Quinoa;

ErrCode
HydroException::handleException(Driver* driver)
//******************************************************************************
//  Handle HydroException
//! \author J. Bakosi
//******************************************************************************
{
  // Contribute to error message
  m_message = HydroMsg[static_cast<int>(m_except)] + m_message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
