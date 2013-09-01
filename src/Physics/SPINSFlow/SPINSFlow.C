//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.C
  \author    J. Bakosi
  \date      Sun 01 Sep 2013 02:27:22 PM MDT
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
                     QuinoaControl* const control,
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
SPINSFlow::init()
//******************************************************************************
//  Initialize the physics
//! \author  J. Bakosi
//******************************************************************************
{
}
