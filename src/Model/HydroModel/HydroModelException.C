//******************************************************************************
/*!
  \file      src/Model/HydroModel/HydroModelException.C
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 10:38:29 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model exception handler
  \details   Hydro model exception handler
*/
//******************************************************************************

#include <HydroModelException.h>

using namespace Quinoa;

ErrCode
HydroModelException::handleException(Driver* driver)
//******************************************************************************
//  Handle HydroModelException
//! \author J. Bakosi
//******************************************************************************
{
  // Contribute to error message
  m_message = HydroModelMsg[static_cast<int>(m_except)] + m_message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
