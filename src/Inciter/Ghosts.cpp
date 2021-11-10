// *****************************************************************************
/*!
  \file      src/Inciter/Ghosts.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Definitions file for generating ghost data structures
  \details   Definitions file for asynchronous distributed
             ghost data structures using Charm++.
*/
// *****************************************************************************

#include "Ghosts.hpp"
#include "DerivedData.hpp"
#include "Vector.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::Ghosts;

Ghosts::Ghosts()
//:
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
}

#include "NoWarning/ghosts.def.h"
