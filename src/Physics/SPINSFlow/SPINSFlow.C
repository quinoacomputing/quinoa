//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.C
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 11:11:30 AM MDT
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

SPINSFlow::SPINSFlow(Memory* const memory,
                     Paradigm* const paradigm,
                     const QuinoaControl& control,
                     Timer* const timer,
                     const QuinoaPrint& print,
                     const std::string& filename) :
  Physics(memory, paradigm, control, timer, print),
  m_filename(filename)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object
//! \param[in]  timer    Timer object pointer
//! \param[in]  filename Mesh filename
//! \author  J. Bakosi
//******************************************************************************
{
}

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
