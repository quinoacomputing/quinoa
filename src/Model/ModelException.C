//******************************************************************************
/*!
  \file      src/Model/ModelException.C
  \author    J. Bakosi
  \date      Fri Apr 26 12:57:24 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     ModelException
  \details   ModelException
*/
//******************************************************************************

#include <ModelException.h>

using namespace Quinoa;

ErrCode
ModelException::handleException(Driver* const driver)
//******************************************************************************
//  Handle ModelException
//! \author J. Bakosi
//******************************************************************************
{
  // Contribute to error message
  m_message = ModelMsg[static_cast<int>(m_except)] + m_message;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
