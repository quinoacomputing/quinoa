//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.C
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 05:02:24 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Standalone-Particle Incompressible Navier-Stokes Flow
  \details   Standalone-Particle Incompressible Navier-Stokes Flow
*/
//******************************************************************************

#include <sstream>

#include <Memory.h>
#include <QuinoaControl.h>
#include <SPINSFlow.h>

using namespace quinoa;

void
SPINSFlow::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
}

void
SPINSFlow::init()
//******************************************************************************
//  Initialize the physics
//! \author  J. Bakosi
//******************************************************************************
{
}

void
SPINSFlow::echo()
//******************************************************************************
//  Echo information on standalone-particle incompressible Navier-Stokes physics
//! \author J. Bakosi
//******************************************************************************
{
  //! Echo information on physics in general
  Physics::echo();

  //! Echo information on standalone-particle Navier-Stokes physics
}
