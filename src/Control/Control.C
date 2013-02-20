//******************************************************************************
/*!
  \file      src/Control/Control.C
  \author    J. Bakosi
  \date      Tue 19 Feb 2013 07:15:29 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Main control category
  \details   Main control category
*/
//******************************************************************************

#include <limits>
#include <tuple>

#include <Control.h>

using namespace Quinoa;
using namespace control;

Control::Control() : m_data(DEFAULTS)
//******************************************************************************
//  Constructor
//  \details Initialize defaults of everything held in Control, some of these
//           will be populated by the parser depending on the control file.
//! \author  J. Bakosi
//******************************************************************************
{
}
