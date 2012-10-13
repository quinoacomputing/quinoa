//******************************************************************************
/*!
  \file      src/Base/RandomException.C
  \author    J. Bakosi
  \date      Thu 11 Oct 2012 08:47:09 PM EDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RandomException class definition
  \details   RandomException class definition
*/
//******************************************************************************

#include <iostream>

#include <RandomException.h>

using namespace Quinoa;

ErrCode
RandomException::handleException(Driver* driver)
//******************************************************************************
//  Handle RandomException
//! \author J. Bakosi
//******************************************************************************
{
  // Output message
  cerr << "Random number generator exception: "
       << RndMsg[static_cast<Int>(m_except)]
       << endl;

  // Handle Exception (criticality)
  return Exception::handleException(driver);
}
