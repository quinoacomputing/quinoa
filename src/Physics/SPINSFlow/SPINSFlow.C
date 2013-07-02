//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.C
  \author    J. Bakosi
  \date      Tue Jul  2 16:06:22 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Standalone-Particle Incompressible Navier-Stokes Flow
  \details   Standalone-Particle Incompressible Navier-Stokes Flow
*/
//******************************************************************************

#include <sstream>

#include <Memory.h>
#include <Control.h>
#include <SPINSFlow.h>

using namespace Quinoa;

SPINSFlow::SPINSFlow(Memory* const memory,
                     Paradigm* const paradigm,
                     Control* const control,
                     Timer* const timer,
                     const std::string& filename) :
  Physics(memory, paradigm, control, timer),
  m_filename(filename)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
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
SPINSFlow::echo() const
//******************************************************************************
//  Echo information on the physics
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
