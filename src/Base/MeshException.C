//******************************************************************************
/*!
  \file      src/Base/MeshException.C
  \author    J. Bakosi
  \date      Wed 12 Sep 2012 06:26:46 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshException class definition
  \details   MeshException class definition
*/
//******************************************************************************

#include <iostream>

#include <MeshException.h>

using namespace Quinoa;

ErrorCode
MeshException::handleException(Driver* driver)
//******************************************************************************
//  Handle MeshException
//! \author J. Bakosi
//******************************************************************************
{
  // Output message
  cerr << "Mesh exception: " << MeshMessage[m_exception];
  if (m_filename.size()) cerr << m_filename;
  cerr << endl;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
