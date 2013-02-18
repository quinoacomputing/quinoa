//******************************************************************************
/*!
  \file      src/Control/Control.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 03:17:00 PM MST
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

Control::Control() :      // defaults
  m_data( std::make_tuple(static_cast<string>(""),     // title
                          PhysicsType::NO_PHYSICS,     // physics
                          HydroType::NO_HYDRO,         // hydro
                          MixType::NO_MIX,             // mix
                          numeric_limits<int>::max(),  // nstep
                          0.0,                         // term
                          1.0,                         // dt
                          1,                           // nscalar
                          1,                           // npar
                          1))                          // echo
//******************************************************************************
//  Constructor
//  \details Initialize defaults of everything held in Control, some of these
//           will be populated by the parser depending on the control file.
//! \author  J. Bakosi
//******************************************************************************
{
}
