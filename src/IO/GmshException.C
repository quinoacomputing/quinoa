//******************************************************************************
/*!
  \file      src/Base/GmshException.C
  \author    J. Bakosi
  \date      Fri 07 Sep 2012 01:17:58 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     GmshException class definition
  \details   GmshException class definition
*/
//******************************************************************************

#include <iostream>

#include <GmshException.h>

using namespace Quinoa;

ErrorCode
GmshException::handleException(Driver* driver)
//******************************************************************************
//  Handle GmshException
//! \author J. Bakosi
//******************************************************************************
{
  // Output message
  cerr << "Gmsh exception: " << GmshMessage[m_exception] << endl;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
