//******************************************************************************
/*!
  \file      src/Base/IOException.C
  \author    J. Bakosi
  \date      Fri Sep 14 17:31:38 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     IOException class definition
  \details   IOException class definition
*/
//******************************************************************************

#include <iostream>

#include <IOException.h>

using namespace Quinoa;

ErrCode
IOException::handleException(Driver* driver)
//******************************************************************************
//  Handle IOException
//! \author J. Bakosi
//******************************************************************************
{
  // Output message
  cerr << "IO exception: " << IOMsg[static_cast<Int>(m_except)];
  if (m_filename.size()) cerr << m_filename;
  cerr << endl;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
