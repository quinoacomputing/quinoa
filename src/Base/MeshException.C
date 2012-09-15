//******************************************************************************
/*!
  \file      src/Base/MeshException.C
  \author    J. Bakosi
  \date      Fri Sep 14 17:36:41 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshException class definition
  \details   MeshException class definition
*/
//******************************************************************************

#include <iostream>

#include <MeshException.h>

using namespace Quinoa;

ErrCode
MeshException::handleException(Driver* driver)
//******************************************************************************
//  Handle MeshException
//! \author J. Bakosi
//******************************************************************************
{
  // Output message
  cerr << "Mesh exception: " << MeshMsg[static_cast<Int>(m_except)];
  if (m_filename.size()) cerr << m_filename;
  cerr << endl;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
