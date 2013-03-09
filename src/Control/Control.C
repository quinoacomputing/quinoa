//******************************************************************************
/*!
  \file      src/Control/Control.C
  \author    J. Bakosi
  \date      Sat 09 Mar 2013 11:38:57 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Main control category
  \details   Main control category
*/
//******************************************************************************

#include <limits>
#include <tuple>

#include <Control.h>

using namespace Quinoa;

Control::Control() : m_data(control::DEFAULTS)
//******************************************************************************
//  Constructor
//  \details Initialize defaults of everything held in Control, some of these
//           will be populated by the parser depending on the control file.
//! \author  J. Bakosi
//******************************************************************************
{
}
