// *****************************************************************************
/*!
  \file      src/Main/LBSwitch.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare group for switching on/off load balancing
  \details   Charm++ chare group for switching on/off load balancing.
*/
// *****************************************************************************

#include <iostream>

#include "LBSwitch.h"
#include "Print.h"

using tk::LBSwitch;

LBSwitch::LBSwitch( bool verbose )
// *****************************************************************************
//  Constructor: turn on automatic load balancing
//! \param[in] verbose True if user selected verbose mode
// *****************************************************************************
{
  TurnManualLBOff();

  if (CkMyPe() == 0)
    Print( verbose ? std::cout : std::clog ).diag( "Load balancing on" );
}

void
LBSwitch::off()
// *****************************************************************************
//  Turn off automatic load balancing
//! \details Since this is a [procinit] routine, the runtime system executes the
//!   routine exactly once on every PE early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  TurnManualLBOn();

  if (CkMyPe() == 0)
    Print( std::cout ).diag( "Load balancing off" );
}

#include "NoWarning/lbswitch.def.h"
