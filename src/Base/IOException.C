//******************************************************************************
/*!
  \file      src/Base/IOException.C
  \author    J. Bakosi
  \date      Thu 13 Sep 2012 04:29:05 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     IOException class definition
  \details   IOException class definition
*/
//******************************************************************************

#include <iostream>

#include <IOException.h>

using namespace Quinoa;

ErrorCode
IOException::handleException(Driver* driver)
//******************************************************************************
//  Handle IOException
//! \author J. Bakosi
//******************************************************************************
{
  // Output message
  cerr << "IO exception: " << IOMessage[m_exception];
  if (m_filename.size()) cerr << m_filename;
  cerr << endl;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
