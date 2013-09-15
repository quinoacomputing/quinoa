//******************************************************************************
/*!
  \file      src/Main/Driver.C
  \author    J. Bakosi
  \date      Sat 14 Sep 2013 08:02:07 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class definition
  \details   Driver base class definition
*/
//******************************************************************************

#include <Driver.h>
#include <Timer.h>
#include <Exception.h>

using namespace quinoa;

Driver::Driver()
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // Instantiate timer object
  m_timer = std::unique_ptr<Timer>(new (std::nothrow) Timer);
  ErrChk(m_timer != nullptr, ExceptType::FATAL,
         "Cannot allocate memory for timer object");
}
