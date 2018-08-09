// *****************************************************************************
/*!
  \file      src/Base/ChareState.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare state collector group
  \details   Charm++ chare state collectory group used for debugging.
*/
// *****************************************************************************

#include "ChareState.h"

using tk::ChareState;

ChareState::ChareState()
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
}

void ChareState::collect() const
// *****************************************************************************
//  Collect chare state
// *****************************************************************************
{
}

#include "NoWarning/charestate.def.h"
